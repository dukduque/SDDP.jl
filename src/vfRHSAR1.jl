module vfRHSAR1

struct RHSAR1CutOracle <: AbstractCutOracle
    cuts::Vector{CutAR1}
end
RHSAR1CutOracle() = RHSAR1CutOracle(Cut[])

struct ValueFunctionRHSAR1{C<:RHSAR1CutOracle} <: AbstractValueFunction
    cutmanager::C
    theta::JuMP.Variable
end
ValueFunctionRHSAR1(cutoracle=RHSAR1CutOracle()) = DefaultValueFunction(cutoracle, JuMP.Variable(JuMP.Model(), 0))


function SDDP.forwardpass!(m::SDDPModel{DanielValueFunction}, settings::SDDP.Settings, solutionstore=nothing)
    ... the forward pass algorithm ...
end

end
