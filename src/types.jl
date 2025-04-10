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

function apply_limiter end

function apply_limiter(limiter::VanLeer, r)
    (r + abs(r)) / (1 + abs(r))
end

function apply_limiter(limiter::VanAlbada, r)
    (r^2 + r) / (r^2 + 1)
end

function apply_limiter(limiter::Minmod, r)
    max(0, min(1, r))
end

function apply_limiter(limiter::Superbee, r)
    max(0, min(1, 2*r), min(2, r))
end

function apply_limiter(limiter::Sweby, r)
    max(0, min(limiter.β*r, 1), min(r, limiter.β))
end

function apply_limiter(limiter::UMIST, r)
    max(0, min(2r, (1+3r)/4, (3+r)/4, 2))
end

function apply_limiter(limiter::UD, r)
    0.0
end

function apply_limiter(limiter::QUICKlimiter, r)
    max(0, min(2r, (3+r)/4, 2))
end
