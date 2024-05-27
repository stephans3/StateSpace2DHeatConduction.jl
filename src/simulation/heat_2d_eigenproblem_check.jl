#=
    Checking the eigenvalue problem of analytically computed eigenvalues and eigenvectors of a 2-dimension heat conduction problem. 
=#


Nx = 20;
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

#D2[1+Nx,1:Nx:1+2Nx]
1+(Ny-1)*Nx
for i=1:Nx
    for j=2:Ny-1
        di = (j-1)*Nx
        D2[i+di,i+di-Nx:Nx:(i+di)+Nx] =  [1,-2,1]
        #D2[i+(j-1)*Nx,i-1+(j-1)*Nx : i+1+(j-1)*Nx] = [1,-2,1];
    end
    D2[i,i:Nx:i+Nx] = [-1,1]
    D2[i+(Ny-1)*Nx,i+(Ny-2)*Nx:Nx:i+(Ny-1)*Nx] = [1,-1]
end

L = 0.2;
W = 0.1;

λ₁, λ₂ = 40, 60;
ρ = 8000;
c = 400;

α₁ = λ₁ / (ρ*c);
α₂ = λ₂ / (ρ*c);

Δx = L/Nx;
Δy = W/Ny;

A = α₁*D1/Δx^2 + α₂*D2/Δy^2;

# Eigenvalues
μdata = zeros(Nx, Ny)

for j=1:Nx, m=1:Ny
    μdata[j,m] = -2α₁*(1-cos((j-1)*π / Nx))/Δx^2 -2α₂*(1-cos((m-1)*π / Ny))/Δy^2
end

μ_ana = reshape(μdata, Nx*Ny)

ψdata = zeros(Nx,Ny,Nx,Ny);

for j = 1:Nx, m = 1:Ny
    for n1=1:Nx, n2=1:Ny
        ψdata[j,m,n1,n2] = cospi.((j-1)*(2n1-1)/(2Nx))*cospi.((m-1)*(2n2-1)/(2Ny))
    end
end

ψ_ana = reshape(ψdata, Nx*Ny, Nx*Ny)

using LinearAlgebra
μ_num, ψ_num =  eigen(A)

# Check Eigenvalues
round.(sort(μ_ana),digits=10) .≈ round.(μ_num, digits=10)


# Check Eigenvalue problem of numerical values
round.(A * ψ_num,digits=10) .≈ round.(ψ_num .* μ_num',digits=10) 

# Check Eigenvalue problem of analytical values
round.(A * inv(ψ_ana),digits=10) .≈ round.(inv(ψ_ana) .* μ_ana', digits=10)

# Build A matrix with numerical Eigenvalues + Eigenvectors
A_num = ψ_num * diagm(μ_num) * inv(ψ_num)

# Build A matrix with analytical Eigenvalues + Eigenvectors
A_ana = inv(ψdata1) * diagm(μdata1) * ψdata1

# Check correctness approximately
A .≈ round.(A_num, digits=10)
A .≈ round.(A_ana, digits=10)

