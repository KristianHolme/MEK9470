abstract type AbstractDiscretization end
struct CentralDiff <: AbstractDiscretization end
struct Upwind <: AbstractDiscretization end
struct Hybrid <: AbstractDiscretization end
struct QUICK <: AbstractDiscretization end