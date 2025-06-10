using DrWatson
@quickactivate :MEK9470
# Run all limiters
for limiter in [Sweby(1.5), VanLeer(), VanAlbada(), Minmod(), Superbee(), UMIST(), QUICKlimiter(), UD(), CD()]
    solve_up(limiter, 512, Res=[1, 10, 20, 40, 70, 100, 150, 200, 250, 350, 500, 750, 1000, 1500, 2000, 3000, 4000, 5000])
end