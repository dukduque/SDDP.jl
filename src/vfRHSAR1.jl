module vfRHSAR1

using SDDP


struct DanielValueFunction <: SDDP.AbstractValueFunction
    ... things you need ...
end

function SDDP.forwardpass!(m::SDDPModel{DanielValueFunction}, settings::SDDP.Settings, solutionstore=nothing)
    ... the forward pass algorithm ...
end

end
