using DrWatson
@quickactivate :MEK9470
# Run all limiters
for limiter in [VanLeer(), VanAlbada(), Minmod(), Superbee(), UMIST(), QUICKlimiter(), UD(), CD()]
    solve_up(limiter, 128)
end