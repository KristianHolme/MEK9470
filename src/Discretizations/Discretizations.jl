module Discretizations

export AbstractDiscretization, AbstractLimiter
export CentralDiff, Upwind, Hybrid, QUICK, TVD
export UD, VanLeer, VanAlbada, Minmod, Superbee, Sweby, UMIST, QUICKlimiter
export apply_limiter

abstract type AbstractDiscretization end
abstract type AbstractLimiter end

struct CentralDiff <: AbstractDiscretization end
struct Upwind <: AbstractDiscretization end
struct Hybrid <: AbstractDiscretization end
struct QUICK <: AbstractDiscretization end
struct TVD <: AbstractDiscretization
    limiter::AbstractLimiter
    function TVD(limiter::AbstractLimiter=VanLeer())
        return new(limiter)
    end
end



struct UD <: AbstractLimiter end
struct VanLeer <: AbstractLimiter end
struct VanAlbada <: AbstractLimiter end
struct Minmod <: AbstractLimiter end
struct Superbee <: AbstractLimiter end
struct Sweby <: AbstractLimiter
    β::Float64
    function Sweby(β::Float64=1.5)
        return new(β)
    end
end
struct UMIST <: AbstractLimiter end
struct QUICKlimiter <: AbstractLimiter end

function apply_limiter(limiter::AbstractLimiter, r::Real)
    limiter(r)
end

function (limiter::VanLeer)(r::Real)
    (r + abs(r)) / (1 + abs(r))
end

function (limiter::VanAlbada)(r::Real)
    (r^2 + r) / (r^2 + 1)
end

function (limiter::Minmod)(r::Real)
    max(0, min(1, r))
end

function (limiter::Superbee)(r::Real)
    max(0, min(1, 2 * r), min(2, r))
end

function (limiter::Sweby)(r::Real)
    max(0, min(limiter.β * r, 1), min(r, limiter.β))
end

function (limiter::UMIST)(r::Real)
    max(0, min(2r, (1 + 3r) / 4, (3 + r) / 4, 2))
end

function (limiter::UD)(r::Real)
    zero(r)
end

function (limiter::QUICKlimiter)(r::Real)
    max(0, min(2r, (3 + r) / 4, 2))
end

end