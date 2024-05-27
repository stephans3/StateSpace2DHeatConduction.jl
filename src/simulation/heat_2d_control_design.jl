#=
    Simulation and controller design for a 2-dimensional heat conduction problem.
=#


# Diffusion matrices
Nx = 30;
Ny = 10;

D1 = zeros(Int64,Nx*Ny, Nx*Ny)
D2 = zeros(Int64,Nx*Ny, Nx*Ny)

for j=1:Ny
    di = (j-1)*Nx
    for i=2:Nx-1
        D1[i+di,i-1+di : i+1+di] = [1,-2,1];
    end
    D1[1+di,1+di:2+di] = [-1,1]
    D1[Nx+di,Nx-1+di:Nx+di] = [1,-1]
end

D1

for i=1:Nx
    for j=2:Ny-1
        di = (j-1)*Nx
        D2[i+di,i+di-Nx:Nx:(i+di)+Nx] =  [1,-2,1]
    end
    D2[i,i:Nx:i+Nx] = [-1,1]
    D2[i+(Ny-1)*Nx,i+(Ny-2)*Nx:Nx:i+(Ny-1)*Nx] = [1,-1]
end


# Material data
λ₁, λ₂ = 40, 60;
ρ = 8000;
c = 400;

α₁ = λ₁ / (ρ*c);
α₂ = λ₂ / (ρ*c);

# Spatial Discretization
L = 0.3;
W = 0.1;

Δx = L/Nx;
Δy = W/Ny;

# System matrix A
A = α₁*D1/Δx^2 + α₂*D2/Δy^2;

# Eigenvalues of A
μdata = zeros(Nx, Ny)
for j=1:Nx, m=1:Ny
    μdata[j,m] = -2α₁*(1-cos((j-1)*π / Nx))/Δx^2 -2α₂*(1-cos((m-1)*π / Ny))/Δy^2
end

μdata1 = reshape(μdata, Nx*Ny)

# Eigenvectors of A
ψdata = zeros(Nx,Ny,Nx,Ny);

for j = 1:Nx, m = 1:Ny
    for n1=1:Nx, n2=1:Ny
        ψdata[j,m,n1,n2] = cospi.((j-1)*(2n1-1)/(2Nx))*cospi.((m-1)*(2n2-1)/(2Ny))
    end
end

ψdata1 = reshape(ψdata, Nx*Ny, Nx*Ny)
ψdata1_inv = inv(ψdata1);

# Time-Sampled System:
# x(n+1) = Ad x(n) + Bd u(n)
Ts = 2
using LinearAlgebra
Ad = ψdata1_inv * diagm(exp.(μdata1*Ts)) * ψdata1

function io_approx(x,n,Nio; M=1, ν=1)
    xc = L*(1+2*(n-1)) / (2*Nio)

    if x > (n-1)*L/Nio && x < n*L/Nio
        return  exp(-(M*(x-xc))^(2ν))
    else
        return 0
    end
end

# Grid of x Points
xgrid = Δx/2 : Δx : L-Δx/2

# Sensors
N_sens = 3
c1 = io_approx.(xgrid,1,N_sens; M=20, ν=1)
c2 = io_approx.(xgrid,2,N_sens; M=20, ν=1)
c3 = io_approx.(xgrid,3,N_sens; M=20, ν=1)

C = vcat(zeros(Int64, Nx*(Ny-1),3), hcat(c1/sum(c1),c2/sum(c2),c3/sum(c3)))'

# Actuators
N_act = 3

b1 = io_approx.(xgrid,1,N_act; M=40, ν=2)
b2 = io_approx.(xgrid,2,N_act; M=40, ν=2)
b3 = io_approx.(xgrid,3,N_act; M=40, ν=2)

# Input matrix
B = vcat(hcat(b1, b2, b3) / (Δy * c * ρ), zeros(Int64, Nx*(Ny-1),3))

function int_B(t,n)
    return ψdata1_inv * diagm(exp.(μdata1*(Ts-t))) * ψdata1 * B[:,n]
end
    
using FastGaussQuadrature
x,w = FastGaussQuadrature.gausslegendre(1000)
p1 = Ts/2

# Discrete input matrix
Bd = zeros(Nx*Ny,3)

intB1 = hcat(int_B.(p1*x .+ p1,1)...)'
intB2 = hcat(int_B.(p1*x .+ p1,2)...)'
intB3 = hcat(int_B.(p1*x .+ p1,3)...)'

Bd[:,1] = (p1*w'*intB1)[:]
Bd[:,2] = (p1*w'*intB2)[:]
Bd[:,3] = (p1*w'*intB3)[:]

# Controller Design with
# Linear Quadratic Regulator
Rw = diagm(ones(N_act));
Qw = 1e5*diagm(ones(Nx*Ny));
Sw = zeros(Nx*Ny, N_act);

# Solving Riccati Equation
using MatrixEquations
P,evals_cl,K = ared(Ad, Bd, Rw, Qw, Sw)

# Closed-loop System matrix
Adcl = Ad-Bd*K

# Static feedforward filter
Wfw = -inv(C*inv(Adcl- I)*Bd)


# Simulation of closed-loop system
Tf = 2000;
tgrid = 0 : Ts : Tf;
num_t = length(tgrid)

Θinit = 10ones(Nx,Ny)
Θinit .= 10*sinpi.(2xgrid / L) * ones(1,Ny);

# Temperature data
Θdata = zeros(Nx*Ny,num_t+1)
Θdata[:,1] = reshape(Θinit, Nx*Ny);

input_data = zeros(3,num_t)

# Input function
function input_u(t,x)
    # Reference
    ref = 5*ones(N_sens);
    return -K*x + Wfw*ref
end

# θ(n+1) = Ad θ(n) +  Bd u(n)
# u(n) = -K θ(n) + W r(n)
for (it, t) in enumerate(tgrid)
    Θdata[:, it+1] = Ad * Θdata[:,it] + B*input_u(0,Θdata[:,it])
    input_data[:,it] = input_u(0,Θdata[:,it])
end



Θdata
input_data

y_out = (C*Θdata)'

using CairoMakie
begin
    fig = Figure(size=(800,600),fontsize=20)
    ax = Axis3(fig[1,1], azimuth = 5pi/4, 
                xlabel = "Time t in [s]", ylabel = "Position x in [m]", zlabel = "Temperature", 
                xlabelsize = 24, ylabelsize = 24, zlabelsize = 24,)

    #surface!(ax, tgrid[1:100:end], xgrid, Θdata[(Ny-1)*Nx+1:Ny*Nx,1:100:end], colormap = :plasma)            
    surface!(ax, tgrid[1:10:end], xgrid, Θdata[(Ny-1)*Nx+1:end,1:10:end-1]', colormap = :plasma)            
    fig
    save("results/figures/"*"temp_north_3d.pdf", fig,pt_per_unit = 1)    
end




begin
    fig1 = Figure(size=(800,600),fontsize=20)
    ax1 = Axis(fig1[1, 1], xlabel = "Time t in [s]", ylabel = L"Input Signal $\times 10^{3}$", ylabelsize = 24,
        xlabelsize = 24, xgridstyle = :dash, ygridstyle = :dash, 
        xtickalign = 1., xticksize = 10, 
        xminorgridvisible = true, xminorticksvisible = true, xminortickalign = 1,
        yminorgridvisible = true, yminorticksvisible = true, yminortickalign = 1,
        ytickalign = 1, yticksize = 10, xlabelpadding = 0)
    
    input_data_scaled = input_data * 1e-3;
    ax1.xticks = 0 : 200 : tgrid[end];    
    ax1.yticks = -5 : 5 : 35;
    lines!(tgrid, input_data_scaled[1,:]; linestyle = :dot,     linewidth = 4, label = L"Input 1$")
    lines!(tgrid, input_data_scaled[2,:]; linestyle = :solid,    linewidth = 4, label = L"Input 2$")
    lines!(tgrid, input_data_scaled[3,:]; linestyle = :dash, linewidth = 4, label = L"Input 3$")
    axislegend(; position = :rt, backgroundcolor = (:grey90, 0.1));

    fig1
    save("results/figures/"*"input_signal.pdf", fig1, pt_per_unit = 1)    
  end


begin
    fig1 = Figure(size=(800,600),fontsize=20)
    ax1 = Axis(fig1[1, 1], xlabel = "Time t in [s]", ylabel = "Temperature in [K]", ylabelsize = 24,
        xlabelsize = 24, xgridstyle = :dash, ygridstyle = :dash, 
        xtickalign = 1., xticksize = 10, 
        xminorgridvisible = true, xminorticksvisible = true, xminortickalign = 1,
        yminorgridvisible = true, yminorticksvisible = true, yminortickalign = 1,
        ytickalign = 1, yticksize = 10, xlabelpadding = 0)
    
    ax1.xticks = 0 : 200 : Tf;    
    ax1.yticks = -7.5 : 2.5 : 7.5;
    lines!(tgrid, y_out[1:end-1,1]; linestyle = :dot,     linewidth = 4, label = "Sensor 1")
    lines!(tgrid, y_out[1:end-1,2]; linestyle = :solid,    linewidth = 4, label = "Sensor 2")
    lines!(tgrid, y_out[1:end-1,3]; linestyle = :dash, linewidth = 4, label = "Sensor 3")
    #scatter!(t_ref, ref_data[1:50:end]; markersize = 15, marker = :diamond, color=:black, label = "Reference")
    axislegend(; position = :rt, backgroundcolor = (:grey90, 0.1));

    fig1
  
    save("results/figures/"*"temp_sensing.pdf", fig1, pt_per_unit = 1)     
end









using Plots
plot((C*tempdata')', legend=false)
plot(input_data, legend=false)

heatmap(reshape(tempdata[1, :], Nx, Ny))

heatmap(reshape(tempdata[200, :], Nx, Ny))


for it = 1:20:2001#length(tgrid) # (it, t) in enumerate(tgrid)
    display(heatmap(reshape(tempdata[it, :], Nx, Ny)))
end


