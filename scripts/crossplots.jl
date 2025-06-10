using DrWatson
@quickactivate :MEK9470
using CairoMakie
using Dates
using JLD2
##
"""
:velocity_magnitude
:u
:v
:pressure
:vorticity
:streamfunction
"""
N = 256
Re = 4000 #250, 4000
quantity = :streamfunction
limiters = [Sweby(1.5), VanLeer(), VanAlbada(), Minmod(), Superbee(), UMIST(), QUICKlimiter(), UD(), CD()]
fig_crossplot = create_crossplot(N, Re, limiters, quantity=quantity, contours=true)
display(fig_crossplot)
name = "all_limiters"
safesave(plotsdir("crossplots", string(quantity), "N=$(N)", "Re=$(Re)", "$(name).png"), fig_crossplot)
##
limiters = [Sweby(1.5), VanLeer(), Minmod(), Superbee(), UMIST(), QUICKlimiter()]
fig_crossplot = create_crossplot(N, Re, limiters, quantity=quantity)
display(fig_crossplot)
name = "no_Va_UD_CD"
safesave(plotsdir("crossplots", string(quantity), "N=$(N)", "Re=$(Re)", "$(name).png"), fig_crossplot)
##
limiters = [Sweby(1.5), VanAlbada(), UD(), CD()]
fig_crossplot = create_crossplot(N, Re, limiters, quantity=quantity)
display(fig_crossplot)
name = "sweby_vanalbada_UD_CD"
safesave(plotsdir("crossplots", string(quantity), "N=$(N)", "Re=$(Re)", "$(name).png"), fig_crossplot)
##
limiters = [VanAlbada(), UD(), CD()]
fig_crossplot = create_crossplot(N, Re, limiters, quantity=quantity)
display(fig_crossplot)
name = "vanalbada_UD_CD"
safesave(plotsdir("crossplots", string(quantity), "N=$(N)", "Re=$(Re)", "$(name).png"), fig_crossplot)


##
shenfun_refs = Dict(250 => (0.61032, 0.663215))
fig_minima = limiter_minima_comparison_plot(N, Re, limiters,
    # limits=((0.6, 0.63), (0.655, 0.675)),
    limits=((0.51, 0.55), (0.529, 0.564)),
    # limits=(nothing, nothing),
    offset_limits=(0.002, 0.002),
    reference_point=get(shenfun_refs, Re, nothing),
    cell_ticks=true)

##
name = "all_limiters_zoom"
safesave(plotsdir("streamfunction_minima", "N=$(N)", "Re=$(Re)", "$(name).png"), fig_minima)
