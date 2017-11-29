#  Copyright 2017,  Oscar Dowson
#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at http://mozilla.org/MPL/2.0/.
#############################################################################

#=
This problem is a version of the hydro-thermal scheduling problem. The goal is
to operate two hydro-dams in a valley chain over time in the face of inflow
and price uncertainty.

Turbine response curves are modelled by piecewise linear functions which map the
flow rate into a power. These can be controlled by specifying the breakpoints
in the piecewise linear function as the knots in the Turbine struct.

The model can be created using the hydrovalleymodel function. It has a few
keyword arguments to allow automated testing of the library.
`hasstagewiseinflows` determines if the RHS noise constraint should be added.
`hasmarkovprice` determines if the price uncertainty (modelled by a markov
chain) should be added.

In the third stage, the markov chain has some unreachable states to test
some code-paths in the library.

We can also set the sense to :Min or :Max (the objective and bound are
flipped appropriately).
=#
using SDDP, JuMP, Clp

immutable Turbine
    flowknots::Vector{Float64}
    powerknots::Vector{Float64}
end

immutable Reservoir
    min::Float64
    max::Float64
    initial::Float64
    turbine::Turbine
    spill_cost::Float64
    normNoise::Vector{Float64}
end

function hydrovalleymodel(;
        riskmeasure=Expectation(),
        cutoracle=DefaultCutOracle(),
        hasstagewiseinflows::Bool=true,
        hasmarkovprice::Bool=true,
        sense::Symbol=:Max
    )
    T = 3
    valley_chain = [
        Reservoir(0, 200, 200, Turbine([50, 60, 70], [55, 65, 70]), 1000, [-12.626063, -3.704914,  5.216234]),
        Reservoir(0, 200, 200, Turbine([50, 60, 70], [55, 65, 70]), 1000, [10.734489, 22.298040, 33.861590]),
        Reservoir(0, 200, 200, Turbine([50, 60, 70], [55, 65, 70]), 1000, [-11.599395 , 5.652446,22.904287]),
        Reservoir(0, 200, 200, Turbine([50, 60, 70], [55, 65, 70]), 1000, [-8.867397 , 8.238426, 25.344249])]

    turbine(i) = valley_chain[i].turbine

    autoreg = [0.894155785206891 0.187191166096471 -0.128736587791406 0.128227942329335;
                -0.328397404318301 0.695183935103107 -0.0913626214178796 0.0289323912343019;
                0.585025574269405 0.531110268093871 -0.113118385105778 -0.16903444953985;
                0.510092039498543 0.203814073453963 0.0196469487416793 0.150465728201917]


    # Prices[stage, markov state]
    prices = [
        1 2 0;
        2 1 0;
        3 4 0;
        4 2 0;
        3 1 0;
        4 1 0;
        3 1 0;
        4 1 0;
    ]

    # Transition matrix
    if hasmarkovprice
        transition = Array{Float64, 2}[
            [ 1.0 ]',
            [ 0.6 0.4 ],
            [ 0.6 0.4 0.0; 0.3 0.7 0.0]
        ]
    else
        transition = [ones(Float64, (1,1)) for t in 1:T]
    end

    flipobj = (sense == :Max)?1.0:-1.0

    N = length(valley_chain)

    # Initialise SDDP Model
    m = SDDPModel(
                sense           = sense,
                stages          = T,
                objective_bound = flipobj * 1e6,
                markov_transition = transition,
                risk_measure    = riskmeasure,
                cut_oracle      = cutoracle,
                solver          = ClpSolver()
                                        ) do sp, stage, markov_state

        # ------------------------------------------------------------------
        #   SDDP State Variables
        # Level of upper reservoir
        @states(sp, begin
                valley_chain[r].min <= reservoir[r=1:N] <= valley_chain[r].max, reservoir0==valley_chain[r].initial
                s_inflow[r=1:N], s_inflow0==0
            end)
        # ------------------------------------------------------------------
        #   Additional variables

        @variables(sp, begin
            outflow[r=1:N]      >= 0
            spill[r=1:N]        >= 0
            pour[r=1:N]        >= 0
            generation_quantity >= 0 # Total quantity of water
            # Proportion of levels to dispatch on
            0 <= dispatch[r=1:N, level=1:length(turbine(r).flowknots)] <= 1
        end)

        # ------------------------------------------------------------------
        # Constraints
        @constraints(sp, begin
            # flow from upper reservoir
            reservoir[1] == reservoir0[1] + s_inflow[1] - outflow[1] - spill[1] + pour[1]

            # other flows
            flow[i=2:N], reservoir[i] == reservoir0[i] + s_inflow[i] - outflow[i] - spill[i] + pour[i] + outflow[i-1] + spill[i-1]

            # Total quantity generated
            generation_quantity == sum(turbine(r).powerknots[level] * dispatch[r,level] for r in 1:N for level in 1:length(turbine(r).powerknots))

            # ------------------------------------------------------------------
            # Flow out
            turbineflow[r=1:N], outflow[r] == sum(turbine(r).flowknots[level] * dispatch[r, level] for level in 1:length(turbine(r).flowknots))

            # Dispatch combination of levels
            dispatched[r=1:N], sum(dispatch[r, level] for level in 1:length(turbine(r).flowknots)) <= 1
        end)

        # rainfall noises
        for i in 1:N

            if hasstagewiseinflows && stage > 1 # in future stages random inflows
                #Inflow autoreg model
                @rhsnoise(sp, rainfall = valley_chain[i].normNoise, s_inflow[i]-sum(autoreg[i,j]*s_inflow0[j] for j in 1:N)==rainfall)
                #@rhsnoise(sp, rainfall = valley_chain[i].inflows, inflow[i] <= rainfall)
            else # in the first stage deterministic inflow
                #@constraint(sp, s_inflow[i]=autoreg*s_inflow0[i])
                @constraint(sp, s_inflow[i] == 0)
            end
        end

        # ------------------------------------------------------------------
        #   Objective Function
        if hasmarkovprice
            @stageobjective(sp, flipobj * (prices[stage, markov_state]*generation_quantity - sum(valley_chain[i].spill_cost*(spill[i]+pour[i])for i in 1:N)))
        else
            @stageobjective(sp, flipobj * (prices[stage, 1]*generation_quantity - sum(valley_chain[i].spill_cost * spill[i] for i in 1:N)))
        end

    end
end
