using Pkg

Pkg.add("CairoMakie")
Pkg.add("HCubature")
Pkg.add("SphericalHarmonics")
Pkg.add("LegendrePolynomials")
Pkg.add("DifferentialEquations")
Pkg.add("LaTeXStrings")


using Printf
using LaTeXStrings
using CairoMakie
using HCubature
using SphericalHarmonics
using LegendrePolynomials
using DifferentialEquations


function flattenSHArray(shArray, lₘₐₓ)
    l = 0
    m = 0
    coeffVector::Array{ComplexF64} = []
    for a in 0:lₘₐₓ * (lₘₐₓ + 2)
        push!(coeffVector, shArray[(l, m)])
        if m == l
            l+=1
            m = -l
        else
            m+=1
        end
    end
    coeffVector
end

function flattenSHArray(shArray)
    l = 0
    m = 0
    coeffVector::Array{ComplexF64} = []
    for a in 1:size(shArray)[1]
        push!(coeffVector, shArray[(l, m)])
        if m == l
            l+=1
            m = -l
        else
            m+=1
        end
    end
    coeffVector
end

#Funciton to expand a function in spherical harmonics
#This just integrates to work out coefficients
expandSphericalHarmonics2(f, lₘₐₓ) = hcubature(x -> sin(x[1]) * f(x[1], x[2]) .* conj(flattenSHArray(computeYlm(x[1], x[2], lmax = lₘₐₓ))), [0, 0], [π, 2 * π])

function LaplacianSphericalCoefficients(coeffs, lₘₐₓ)
    for l in 0:lₘₐₓ
        for m in -l:l
            coeffs[(l, m)] *= - l * (l + 1)
        end
    end
    coeffs
end

#Calculate the Laplacian of a representation of a function
#in spherical harmonics given the coefficients
function LaplacianSphericalCoefficientsArays(coeffs, lₘₐₓ)
    l = 0
    m = 0
    coeffVector::Array{ComplexF64} = []
    for a in 0:lₘₐₓ * (lₘₐₓ + 2)
        push(coeffVector, - l * (l + 1) *coeffs[a + 1])
        if m == l
            l+=1
            m = -l
        else
            m+=1
        end
    end
    coeffVector
end

function LaplacianSphericalCoefficientsArrays(coeffs)
    l = 0
    m = 0
    coeffVector::Array{ComplexF64} = []
    for a in coeffs
        push!(coeffVector, - l * (l + 1) * a)
        if m == l
            l+=1
            m = -l
        else
            m+=1
        end
    end
    coeffVector
end


#Calculate a function from spherical harmonic represntation coefficients
function functionFromSphericalCoefficients(coeffs)
    l = 0
    m = 0
    lₘₐₓ = 0
    for a in 1:size(coeffs)[1]
        lₘₐₓ = l
        if m == l
            l+=1
            m = -l
        else
            m+=1
        end
    end
    function f(θ, ϕ)
      harmonics = flattenSHArray(computeYlm(θ, ϕ, lmax = lₘₐₓ))
      harmonics' * coeffs
    end
end

function waveSystem(U)
    half = Int(size(U)[1] // 2)
    ϕ = U[begin:half]
    ψ = U[half + 1:end]    
    vcat(ψ, LaplacianSphericalCoefficientsArrays(ϕ))
end

function solveWaveEquation(ϕ₀, ψ₀, N, t₀, t₁)
  ϕ₀SphericalCoeffs = expandSphericalHarmonics2(ϕ₀, N)[1];
  ψ₀SphericalCoeffs = expandSphericalHarmonics2(ψ₀, N)[1];

  alg = DP5()
  U0 = vcat( ϕ₀SphericalCoeffs, ψ₀SphericalCoeffs)
  prob = ODEProblem((U,p,t) -> waveSystem(U), U0, (t₀, t₁))
  solve(prob, alg)
end

for lmax in [2, 4, 6, 8, 10]
  ϕ₀(θ, ϕ) = exp( - θ^2 / 0.4);
  ψ₀(θ, ϕ) = 0;
  sol = solveWaveEquation(ϕ₀, ψ₀, lmax, 0, 10)

  #Only take the first items
  sizA = Int(size(sol.u[1])[1] / 2)

  θs = LinRange(0, π, 100)
  ϕs = LinRange(0,  2 * π, 100)

  xs = [sin(θ) * cos(ϕ) for θ in θs, ϕ in ϕs]

  ys = [sin(θ) * sin(ϕ) for θ in θs, ϕ in ϕs]

  zs = [cos(θ) for θ in θs, ϕ in ϕs]

  firstprofileWhole = sol[1][begin:sizA]
  myFuncFirst = functionFromSphericalCoefficients(firstprofileWhole)

  cols = [real(myFuncFirst(θ, ϕ)) for θ in θs, ϕ in ϕs]

  points = Observable(zs)
  HB_vn = Observable(cols)

  fps = 10
  nframes = size(sol.t)[1]

  scene = surface(xs, ys, points, color = HB_vn, axis=(type=Axis3,), colorrange=(-0.5,1), colormap = reverse(cgrad(:thermal)))

  record(scene, "GaussianWaveSurfaceOnSphereFromScript"*string(lmax)*".gif"; framerate = fps) do io
      for it = 1:nframes
          firstprofileWhole = sol[it][begin:sizA]
          myFuncFirst = functionFromSphericalCoefficients(firstprofileWhole)    
          cols = [real(myFuncFirst(θ, ϕ)) for θ in θs, ϕ in ϕs]
          points[] = zs
          HB_vn[] = cols
          recordframe!(io)
      end
  end

  firstprofileWhole = sol[1][begin:sizA]
  myFuncFirst = functionFromSphericalCoefficients(firstprofileWhole)

  zs = [real(myFuncFirst(θ, ϕ)) for θ in θs, ϕ in ϕs]

  HB_vn = points = Observable(zs)

  fps = 10
  nframes = size(sol)[2]

  scene = surface(θs, ϕs, HB_vn, axis=(type=Axis3,), colorrange=(-1,1), colormap = reverse(cgrad(:thermal)))


  record(scene, "GaussianWaveSurfaceFromScript"*string(lmax)*".gif"; framerate = fps) do io
      for it = 1:nframes
          firstprofileWhole = sol[it][begin:sizA]
          myFuncFirst = functionFromSphericalCoefficients(firstprofileWhole)    
          zs = [real(myFuncFirst(θ, ϕ)) for θ in θs, ϕ in ϕs]
      
          HB_vn[] = zs
          recordframe!(io)
      end
  end

  nslices = 5

  f = Figure(size=(520, nslices * 320))#Figure(size=(400, trunc(Int, 400 * 5)))
  for (k, i) in enumerate(1:trunc(Int, nframes / nslices):nframes)#trunc(Int, nframes / 5 + 1):nframes
    Axis3(f[k, 1], viewmode = :fit)
    firstprofileWhole = sol[i][begin:sizA]
    myFuncFirst = functionFromSphericalCoefficients(firstprofileWhole)    
    zs = [real(myFuncFirst(θ, ϕ)) for θ in θs, ϕ in ϕs]
    surface!(θs, ϕs, zs, colorrange=(-1,1), colormap = reverse(cgrad(:thermal)))
  end

  save("GaussianWaveSurfaceFromScriptTimes"*string(lmax)*".png", f)

  f = Figure(size=(520, nslices * 320))#Figure(size=(400, trunc(Int, 400 * 5)))
  for (k, i) in enumerate(1:trunc(Int, nframes / nslices):nframes)#trunc(Int, nframes / 5 + 1):nframes
    solTime = @sprintf "%.2f" sol.t[i]
#    title=L"Solution at time $t = %$(solTime)$"
    ax = Axis(f[k, 1], title=L"Solution at time $t = %$(solTime)$", xlabel = L"\theta", ylabel = L"\phi")#@sprintf "Solution at t = %.2fs" sol.t[i], xlabel = L"\sum_a^b{xy} + \mathscr{L}", ylabel = L"\sum_a^b{xy} + \mathscr{L}")
    firstprofileWhole = sol[i][begin:sizA]
    myFuncFirst = functionFromSphericalCoefficients(firstprofileWhole)    
    zs = (y -> real(myFuncFirst(y, pi / 2))).(θs)
    lines!(θs, zs)
  end
  save("GaussianWaveSurfaceFromScriptTimesSlice"*string(lmax)*".png", f)
end


function heat(U)
  LaplacianSphericalCoefficientsArrays(U)
end

ϕ₀(θ, ϕ) = exp(- θ^2 / 0.4)
U0 = expandSphericalHarmonics2(ϕ₀, 15)[1]
prob = ODEProblem((U,p,t) -> heat(U), U0, (0.0, 10.0))

alg = Rodas5(autodiff = false)
solheat = solve(prob, alg);

θs = LinRange(0, π, 100)
ϕs = LinRange(0,  2 * π, 100)

xs = [sin(θ) * cos(ϕ) for θ in θs, ϕ in ϕs]

ys = [sin(θ) * sin(ϕ) for θ in θs, ϕ in ϕs]

zs = [cos(θ) for θ in θs, ϕ in ϕs]

firstprofileWhole = solheat[1][begin:256]
myFuncFirst = functionFromSphericalCoefficients(firstprofileWhole)

cols = [real(myFuncFirst(θ, ϕ)) for θ in θs, ϕ in ϕs]

points = Observable(zs)
HB_vn = Observable(cols)

fps = 1
nframes = size(solheat.t)[1]

scene = surface(xs, ys, points, color = HB_vn, axis=(type=Axis3,), colorrange=(-0.2,0.9), colormap = reverse(cgrad(:redsblues)))

record(scene, "heatSurfaceOnSphereFromScript.gif"; framerate = fps) do io
  for it = 1:nframes
    firstprofileWhole = solheat[it][begin:256]
    myFuncFirst = functionFromSphericalCoefficients(firstprofileWhole)    
    cols = [real(myFuncFirst(θ, ϕ)) for θ in θs, ϕ in ϕs]
    points[] = zs
    HB_vn[] = cols
    recordframe!(io)
  end
end

nslices = 6

f = Figure(size=(520, div(nslices * 320, 2)))
for (k, i) in enumerate(1:trunc(Int, nframes / nslices):nframes)
  solTime = @sprintf "%.2f" solheat.t[i]
  Axis3(f[div(k - 1, 2) + 1 , (k - 1) % 2 + 1], title=L"Solution at time $t = %$(solTime)$", viewmode = :fit)
  profile = solheat[i]
  solFunc = functionFromSphericalCoefficients(profile)    
  cols = [real(solFunc(θ, ϕ)) for θ in θs, ϕ in ϕs]
  surface!(xs, ys, zs, colorrange=(-0.4,0.7), color = cols, colormap = :thermal)
end
save("HeatTimesOnSphere.png", f)

###
####θs = LinRange(0, π, 100)
####ϕs = LinRange(0,  2 * π, 100)
####
####xs = [sin(θ) * cos(ϕ) for θ in θs, ϕ in ϕs]
####
####ys = [sin(θ) * sin(ϕ) for θ in θs, ϕ in ϕs]
####
####zs = [cos(θ) for θ in θs, ϕ in ϕs]
####
####firstprofileWhole = sol[1][begin:256]
####myFuncFirst = functionFromSphericalCoefficients(firstprofileWhole)
####
####cols = [real(myFuncFirst(θ, ϕ)) for θ in θs, ϕ in ϕs]
####
####points = Observable(zs)
####HB_vn = Observable(cols)
####
####fps = 5
####nframes = size(sol)[2]
####
####scene = surface(xs, ys, points, color = HB_vn, axis=(type=Axis3,), colorrange=(-0.5,0.5), colormap = reverse(cgrad(:redsblues)))
####
####record(scene, "waveSurfaceOnSphereFromScript.gif"; framerate = fps) do io
####    for it = 1:nframes
####        firstprofileWhole = sol[it][begin:256]
####        myFuncFirst = functionFromSphericalCoefficients(firstprofileWhole)    
####        cols = [real(myFuncFirst(θ, ϕ)) for θ in θs, ϕ in ϕs]
####        points[] = zs
####        HB_vn[] = cols
####        recordframe!(io)
####        #(...) 
####    end
####end
####
####
####θs = LinRange(0, π, 100)
####ϕs = LinRange(0, 2 * π, 100)
####
####firstprofileWhole = sol[1][begin:256]
####myFuncFirst = functionFromSphericalCoefficients(firstprofileWhole)
####
####zs = [real(myFuncFirst(θ, ϕ)) for θ in θs, ϕ in ϕs]
####
####HB_vn = points = Observable(zs)
####
####fps = 10
####nframes = size(sol)[2]
####
####scene = surface(θs, ϕs, HB_vn, axis=(type=Axis3,), colorrange=(-0.5,0.5), colormap = reverse(cgrad(:redsblues)))
####
####
####record(scene, "waveSurfaceFromScript.gif"; framerate = fps) do io
####    for it = 1:nframes
####        firstprofileWhole = sol[it][begin:256]
####        myFuncFirst = functionFromSphericalCoefficients(firstprofileWhole)    
####        zs = [real(myFuncFirst(θ, ϕ)) for θ in θs, ϕ in ϕs]
####    
####        HB_vn[] = zs
####        recordframe!(io)
####        #(...) 
####    end
####end
###
###ϕ₀(θ, ϕ) = exp(- 4 * θ^4);
###ψ₀(θ, ϕ) = 0;
###
###ϕ₀SphericalCoeffs = expandSphericalHarmonics2(ϕ₀, 4)[1];
###ψ₀SphericalCoeffs = expandSphericalHarmonics2(ψ₀, 4)[1];
###
###U0 = vcat( flattenSHArray(ϕ₀SphericalCoeffs), flattenSHArray(ψ₀SphericalCoeffs))
###prob = ODEProblem((U,p,t) -> waveSystem(U), U0, (0.0, 10.0))
###sol = solve(prob);
###
###θs = LinRange(0, π, 100)
###ϕs = LinRange(0, 2 * π, 100)
###
###sizA = Int(size(sol[1])[1] / 2)
###
###firstprofileWhole = sol[1][1:sizA]
###myFuncFirst = functionFromSphericalCoefficients(firstprofileWhole)
###
###zs = [real(myFuncFirst(θ, ϕ)) for θ in θs, ϕ in ϕs]
###
###HB_vn = points = Observable(zs)
###
###fps = 10
###nframes = size(sol)[2]
###
###scene = surface(θs, ϕs, HB_vn, axis=(type=Axis3,), colorrange=(-0.5,0.5), colormap = reverse(cgrad(:redsblues)))
###
###
###record(scene, "waveSurfaceFromScript4.gif"; framerate = fps) do io
###    for it = 1:nframes
###        firstprofileWhole = sol[it][begin:sizA]
###        myFuncFirst = functionFromSphericalCoefficients(firstprofileWhole)    
###        zs = [real(myFuncFirst(θ, ϕ)) for θ in θs, ϕ in ϕs]
###    
###        HB_vn[] = zs
###        recordframe!(io)
###        #(...) 
###    end
###end
###
####θs = LinRange(0, π, 100)
####ϕs = LinRange(0, 2 * π, 100)
####
####firstprofileWhole = sol[1][begin:256]
####myFuncFirst = functionFromSphericalCoefficients(firstprofileWhole)
####
####zs = [real(myFuncFirst(θ, ϕ)) for θ in θs, ϕ in ϕs]
####
####HB_vn = points = Observable(zs)
####
####fps = 10
####nframes = size(sol)[2]
####
####scene = surface(θs, ϕs, HB_vn, axis=(type=Axis3,), colorrange=(-0.5,0.5), colormap = reverse(cgrad(:redsblues)))
####
####
####record(scene, "waveSurfaceFromScript.gif"; framerate = fps) do io
####    for it = 1:nframes
####        firstprofileWhole = sol[it][begin:256]
####        myFuncFirst = functionFromSphericalCoefficients(firstprofileWhole)    
####        zs = [real(myFuncFirst(θ, ϕ)) for θ in θs, ϕ in ϕs]
####    
####        HB_vn[] = zs
####        recordframe!(io)
####        #(...) 
####    end
####end
###
####
####ϕ₀(θ, ϕ) = exp(- 4 * (ϕ^4 + (θ - pi / 2)^4));
####ψ₀(θ, ϕ) = 0;
####
####ϕ₀SphericalCoeffs = expandSphericalHarmonics2(ϕ₀, 15)[1];
####ψ₀SphericalCoeffs = expandSphericalHarmonics2(ψ₀, 15)[1];
####
####U0 = vcat( flattenSHArray(ϕ₀SphericalCoeffs, 15), flattenSHArray(ψ₀SphericalCoeffs, 15))
####prob = ODEProblem((U,p,t) -> waveSystem(U), U0, (0.0, 10.0))
####waveSol2 = solve(prob);
####
####θs = LinRange(0, π, 100)
####ϕs = LinRange(0,  2 * π, 100)
####
####xs = [sin(θ) * cos(ϕ) for θ in θs, ϕ in ϕs]
####
####ys = [sin(θ) * sin(ϕ) for θ in θs, ϕ in ϕs]
####
####zs = [cos(θ) for θ in θs, ϕ in ϕs]
####
####firstprofileWhole = waveSol2[1][begin:256]
####myFuncFirst = functionFromSphericalCoefficients(firstprofileWhole)
####
####cols = [real(myFuncFirst(θ, ϕ)) for θ in θs, ϕ in ϕs]
####
####points = Observable(zs)
####HB_vn = Observable(cols)
####
####fps = 10
####nframes = size(waveSol2)[2]
####
####scene = surface(xs, ys, points, color = HB_vn, axis=(type=Axis3,), colorrange=(-0.5,0.5))
####
####record(scene, "waveSurfaceOnSphere2FromScript.gif"; framerate = 25) do io
####    for it = 1:nframes
####        firstprofileWhole = waveSol2[it][begin:256]
####        myFuncFirst = functionFromSphericalCoefficients(firstprofileWhole)    
####        cols = [real(myFuncFirst(θ, ϕ)) for θ in θs, ϕ in ϕs]
####        points[] = zs
####        HB_vn[] = cols
####        recordframe!(io)
####        #(...) 
####    end
####end
####
####



#θs = LinRange(0, π, 100)
#ϕs = LinRange(0, 2 * π, 100)
#
##zs = [ϕ₀(θ, ϕ) for θ in θs, ϕ in ϕs]
##
##
##myFig = surface(θs, ϕs, zs, axis=(type=Axis3,), colorrange=(-1,1), colormap = reverse(cgrad(:redsblues)))
##
##save("sin.png", myFig)
## 
#ϕ₀SphericalCoeffs = expandSphericalHarmonics2(ϕ₀, 15)[1];
#
##println(ϕ₀SphericalCoeffs)
#
#MyFuncInBasis = functionFromSphericalCoefficients(ϕ₀SphericalCoeffs)
#
#newzs = [real(MyFuncInBasis(θ, ϕ)) for θ in θs, ϕ in ϕs]
#
#myFig2 = surface(θs, ϕs, newzs, axis=(type=Axis3,), colorrange=(-1,1), colormap = reverse(cgrad(:redsblues)))
#
#save("sinExpanded.png", myFig2)
#
#exit()


#ϕ₀(θ, ϕ) = cos(θ)^2;
#ψ₀(θ, ϕ) = 0;
#
#sol = solveWaveEquation(ϕ₀, ψ₀, 7, 0, 20)
#
#sizA = Int(size(sol.u[1])[1] / 2)
#
#θs = LinRange(0, π, 100)
#ϕs = LinRange(0,  2 * π, 100)
#
#xs = [sin(θ) * cos(ϕ) for θ in θs, ϕ in ϕs]
#
#ys = [sin(θ) * sin(ϕ) for θ in θs, ϕ in ϕs]
#
#zs = [cos(θ) for θ in θs, ϕ in ϕs]
#
#firstprofileWhole = sol[1][begin:sizA]
#myFuncFirst = functionFromSphericalCoefficients(firstprofileWhole)
#
#cols = [real(myFuncFirst(θ, ϕ)) for θ in θs, ϕ in ϕs]
#
#points = Observable(zs)
#HB_vn = Observable(cols)
#
#fps = 10
#nframes = size(sol.t)[1]
#
#scene = surface(xs, ys, points, color = HB_vn, axis=(type=Axis3,), colorrange=(-0,1), colormap = reverse(cgrad(:thermal)))
#
#record(scene, "cosSquaredWaveSurfaceOnSphereFromScript.gif"; framerate = fps) do io
#    for it = 1:nframes
#        firstprofileWhole = sol[it][begin:sizA]
#        myFuncFirst = functionFromSphericalCoefficients(firstprofileWhole)    
#        cols = [real(myFuncFirst(θ, ϕ)) for θ in θs, ϕ in ϕs]
#        points[] = zs
#        HB_vn[] = cols
#        recordframe!(io)
#    end
#end
#
#firstprofileWhole = sol[1][begin:sizA]
#myFuncFirst = functionFromSphericalCoefficients(firstprofileWhole)
#
#zs = [real(myFuncFirst(θ, ϕ)) for θ in θs, ϕ in ϕs]
#
#HB_vn = points = Observable(zs)
#
#fps = 10
#nframes = size(sol)[2]
#
#scene = surface(θs, ϕs, HB_vn, axis=(type=Axis3,), colorrange=(0,1), colormap = reverse(cgrad(:thermal)))
#
#
#record(scene, "cosSquaredWaveSurfaceFromScript.gif"; framerate = fps) do io
#    for it = 1:nframes
#        firstprofileWhole = sol[it][begin:sizA]
#        myFuncFirst = functionFromSphericalCoefficients(firstprofileWhole)    
#        zs = [real(myFuncFirst(θ, ϕ)) for θ in θs, ϕ in ϕs]
#    
#        HB_vn[] = zs
#        recordframe!(io)
#    end
#end

