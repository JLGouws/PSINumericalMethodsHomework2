using Pkg

Pkg.add("WGLMakie")
Pkg.add("HCubature")
Pkg.add("SphericalHarmonics")
Pkg.add("LegendrePolynomials")
Pkg.add("DifferentialEquations")


using WGLMakie
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

expandSphericalHarmonics2(f, lₘₐₓ) = hcubature(x -> sin(x[1]) * f(x[1], x[2]) .* conj(flattenSHArray(computeYlm(x[1], x[2], lmax = lₘₐₓ))), [0, 0], [π, 2 * π])

function LaplacianSphericalCoefficients(coeffs, lₘₐₓ)
    for l in 0:lₘₐₓ
        for m in -l:l
            coeffs[(l, m)] *= - l * (l + 1)
        end
    end
    coeffs
end

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

chunk(arr, n) = [arr[i:min(i + n - 1, end)] for i in 1:n:length(arr)]

function waveSystem(U)
    half = Int(size(U)[1] // 2)
    ϕ = U[begin:half]
    ψ = U[half + 1:end]    
    vcat(ψ, LaplacianSphericalCoefficientsArrays(ϕ))
end

##ϕ₀(θ, ϕ) = exp(- 4 * θ^4);
##ψ₀(θ, ϕ) = 0;
##
##ϕ₀SphericalCoeffs = expandSphericalHarmonics2(ϕ₀, 15)[1];
##ψ₀SphericalCoeffs = expandSphericalHarmonics2(ψ₀, 15)[1];

##U0 = vcat( flattenSHArray(ϕ₀SphericalCoeffs, 15), flattenSHArray(ψ₀SphericalCoeffs, 15))
##prob = ODEProblem((U,p,t) -> waveSystem(U), U0, (0.0, 10.0))
##sol = solve(prob);

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

function solveWaveEquation(ϕ₀, ψ₀, N, t₀, t₁)
  ϕ₀SphericalCoeffs = expandSphericalHarmonics2(ϕ₀, N)[1];
  ψ₀SphericalCoeffs = expandSphericalHarmonics2(ψ₀, N)[1];


  U0 = vcat( ϕ₀SphericalCoeffs, N, ψ₀SphericalCoeffs, N)
  prob = ODEProblem((U,p,t) -> waveSystem(U), U0, (t₀, t₁))
  solve(prob);
end

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
###function heat(U)
###    LaplacianSphericalCoefficientsArrays(U)
###end
###
###ϕ₀(θ, ϕ) = exp(- θ^2)
###U0 = flattenSHArray(expandSphericalHarmonics2(ϕ₀, 15)[1])
###prob = ODEProblem((U,p,t) -> heat(U), U0, (0.0, 10.0))
###
###alg = Rodas5(autodiff = false)
###solheat = solve(prob, alg);
###
###θs = LinRange(0, π, 100)
###ϕs = LinRange(0,  2 * π, 100)
###
###xs = [sin(θ) * cos(ϕ) for θ in θs, ϕ in ϕs]
###
###ys = [sin(θ) * sin(ϕ) for θ in θs, ϕ in ϕs]
###
###zs = [cos(θ) for θ in θs, ϕ in ϕs]
###
###firstprofileWhole = solheat[1][begin:256]
###myFuncFirst = functionFromSphericalCoefficients(firstprofileWhole)
###
###cols = [real(myFuncFirst(θ, ϕ)) for θ in θs, ϕ in ϕs]
###
###points = Observable(zs)
###HB_vn = Observable(cols)
###
###fps = 1
###nframes = size(solheat.t)[1]
###
###scene = surface(xs, ys, points, color = HB_vn, axis=(type=Axis3,), colorrange=(-0.2,0.9), colormap = reverse(cgrad(:redsblues)))
###
###record(scene, "heatSurfaceOnSphereFromScript.gif"; framerate = fps) do io
###    for it = 1:nframes
###        firstprofileWhole = solheat[it][begin:256]
###        myFuncFirst = functionFromSphericalCoefficients(firstprofileWhole)    
###        cols = [real(myFuncFirst(θ, ϕ)) for θ in θs, ϕ in ϕs]
###        points[] = zs
###        HB_vn[] = cols
###        recordframe!(io)
###        #(...) 
###    end
###end



ϕ₀(θ, ϕ) = cos(8 * θ);
ψ₀(θ, ϕ) = 0;

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

ϕ₀(θ, ϕ) = cos(8 * θ);
ψ₀(θ, ϕ) = 0;

sol = solveWaveEquation(ϕ₀, ψ₀, 10, 0, 10)

sizA = Int(size(sol[1])[1] / 2)

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

fps = 1
nframes = size(sol.t)[1]

scene = surface(xs, ys, points, color = HB_vn, axis=(type=Axis3,), colorrange=(-0.5,0.5), colormap = reverse(cgrad(:redsblues)))

record(scene, "cosWaveSurfaceOnSphereFromScript.gif"; framerate = fps) do io
    for it = 1:nframes
        firstprofileWhole = sol[it][begin:sizA]
        myFuncFirst = functionFromSphericalCoefficients(firstprofileWhole)    
        cols = [real(myFuncFirst(θ, ϕ)) for θ in θs, ϕ in ϕs]
        points[] = zs
        HB_vn[] = cols
        recordframe!(io)
    end
end

θs = LinRange(0, π, 100)
ϕs = LinRange(0, 2 * π, 100)

firstprofileWhole = sol[1][begin:sizA]
myFuncFirst = functionFromSphericalCoefficients(firstprofileWhole)

zs = [real(myFuncFirst(θ, ϕ)) for θ in θs, ϕ in ϕs]

HB_vn = points = Observable(zs)

fps = 1
nframes = size(sol)[2]

scene = surface(θs, ϕs, HB_vn, axis=(type=Axis3,), colorrange=(-0.5,0.5), colormap = reverse(cgrad(:redsblues)))


record(scene, "cosWaveSurfaceFromScript.gif"; framerate = fps) do io
    for it = 1:nframes
        firstprofileWhole = sol[it][begin:sizA]
        myFuncFirst = functionFromSphericalCoefficients(firstprofileWhole)    
        zs = [real(myFuncFirst(θ, ϕ)) for θ in θs, ϕ in ϕs]
    
        HB_vn[] = zs
        recordframe!(io)
    end
end
