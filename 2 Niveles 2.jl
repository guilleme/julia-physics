
using DifferentialEquations

atomo0 = [1.0 + im*0.0, 0.0 + im*0.0] 

omega = 2.5pi
omega0 = 3pi 
omega_rabi = pi

delta = omega - omega0 

tspan = (0.0,4.0) 

function ecuacion(datomo,atomo,p,t) 
    datomo[1] = -omega_rabi*im*atomo[2]*cos(omega*t)
    datomo[2] = -omega_rabi*im*atomo[1]*cos(omega*t) + im*omega0*atomo[2]
end 

prob = ODEProblem(ecuacion,atomo0,tspan) 

sol = solve(prob; save_everystep = true, saveat=0.01) 

ce = sol[2,:]
cg = sol[1,:] 

pe = abs.(.*(ce,conj(ce)))
pg = abs.(.*(cg,conj(cg)))

plot(sol.t,pe, label="Pe") 
plot!(sol.t,pg, label = "Pg", xlabel="Tiempo", ylabel="Probabilidad",title="delta = 0, omega = 1 pi") 


# Comienza la ecuación con la RWA: 

atomo0 = [1.0 + im*0.0, 0.0 + im*0.0] 

tspan = (0.0,4.0) 

function ecuacion(datomo,atomo,p,t) 
    datomo[1] = -0.5*omega_rabi*im*atomo[2]
    datomo[2] = -0.5*omega_rabi*im*atomo[1] + im*delta*atomo[2]
end 

prob = ODEProblem(ecuacion,atomo0,tspan) 

sol = solve(prob; save_everystep = true, saveat=0.01) 

ce = sol[2,:]
cg = sol[1,:] 

pe = abs.(.*(ce,conj(ce)))
pg = abs.(.*(cg,conj(cg)))

plot!(sol.t,pe, label="Pe con RWA") 
plot!(sol.t,pg, label = "Pg con RWA", xlabel="Tiempo", ylabel="Probabilidad",title="delta = - 0.5 pi, omega = 2.5 pi")


savefig("RWA_cond05.png")

gr()



omega = pi
gamma = 0.1*pi
delta = 0.3*pi
ee = 0.0 + im*0.0 
ge = 0.0 + im*0.0 
eg = 0.0 + im*0.0 
gg = 1.0 + im*0.0 
rho0 = [gg,eg,ge,ee]
tspan =(0.0,4.0) 


function evolucion(drho,rho,p,t)
 drho[1] = -0.5*im*omega*(rho[2] - rho[3]) + gamma*rho[4] 
 drho[2] = 0.5*im*omega*(rho[4] - rho[1]) - 0.5*gamma*rho[2] + im*delta*rho[2]
 drho[3] = -0.5*im*omega*(rho[4] - rho[1]) - 0.5*gamma*rho[3] - im*delta*rho[3]
 drho[4] = 0.5*im*omega*(rho[2] - rho[3]) - gamma*rho[4] 
end 

prob2 = ODEProblem(evolucion,rho0,tspan) 

sol2 = solve(prob2; save_everystep = true, saveat=0.01) 

peet = abs.(.*(sol2[4,:],conj(sol2[4,:])))
pggt = abs.(.*(sol2[1,:],conj(sol2[1,:])))
pegt = abs.(.*(sol2[2,:],conj(sol2[2,:])))
pget = abs.(.*(sol2[3,:],conj(sol2[3,:])))


plot(sol2.t,peet, label="Pe")
plot!(sol2.t,pggt, label ="Pg",xlabel="tiempo", ylabel="Probabilidad", title="Gamma = 0.1pi, delta = 0.3pi")


plot()
omega = pi
gamma = 0*pi
delta = 0*pi
ee = 0.0 + im*0.0 
ge = 0.0 + im*0.0 
eg = 0.0 + im*0.0 
gg = 1.0 + im*0.0 
rho0 = [gg,eg,ge,ee]
tspan =(0.0,4.0) 


function evolucion(drho,rho,p,t)
 drho[1] = -0.5*im*omega*(rho[2] - rho[3]) + gamma*rho[4] 
 drho[2] = 0.5*im*omega*(rho[4] - rho[1]) - 0.5*gamma*rho[2] + im*delta*rho[2]
 drho[3] = -0.5*im*omega*(rho[4] - rho[1]) - 0.5*gamma*rho[3] - im*delta*rho[3]
 drho[4] = 0.5*im*omega*(rho[2] - rho[3]) - gamma*rho[4] 
end 

prob2 = ODEProblem(evolucion,rho0,tspan) 

sol2 = solve(prob2; save_everystep = true, saveat=0.01) 

peet = abs.(sol2[4,:])

plot!(sol2.t,peet, label ="Pe, Γ = 0, δ = 0",xlabel="tiempo", ylabel="Probabilidad", title="Prob. vs Tiempo")

omega = pi
gamma = 0.2*pi
delta = 0*pi
ee = 0.0 + im*0.0 
ge = 0.0 + im*0.0 
eg = 0.0 + im*0.0 
gg = 1.0 + im*0.0 
rho0 = [gg,eg,ge,ee]
tspan =(0.0,4.0) 


function evolucion(drho,rho,p,t)
 drho[1] = -0.5*im*omega*(rho[2] - rho[3]) + gamma*rho[4] 
 drho[2] = 0.5*im*omega*(rho[4] - rho[1]) - 0.5*gamma*rho[2] + im*delta*rho[2]
 drho[3] = -0.5*im*omega*(rho[4] - rho[1]) - 0.5*gamma*rho[3] - im*delta*rho[3]
 drho[4] = 0.5*im*omega*(rho[2] - rho[3]) - gamma*rho[4] 
end 

prob2 = ODEProblem(evolucion,rho0,tspan) 

sol2 = solve(prob2; save_everystep = true, saveat=0.01) 

peet = abs.(sol2[4,:])

plot!(sol2.t,peet, label ="Pe, Γ = 0.2π, δ = 0π")

omega = pi
gamma = 0.5*pi
delta = 0*pi
ee = 0.0 + im*0.0 
ge = 0.0 + im*0.0 
eg = 0.0 + im*0.0 
gg = 1.0 + im*0.0 
rho0 = [gg,eg,ge,ee]
tspan =(0.0,4.0) 


function evolucion(drho,rho,p,t)
 drho[1] = -0.5*im*omega*(rho[2] - rho[3]) + gamma*rho[4] 
 drho[2] = 0.5*im*omega*(rho[4] - rho[1]) - 0.5*gamma*rho[2] + im*delta*rho[2]
 drho[3] = -0.5*im*omega*(rho[4] - rho[1]) - 0.5*gamma*rho[3] - im*delta*rho[3]
 drho[4] = 0.5*im*omega*(rho[2] - rho[3]) - gamma*rho[4] 
end 

prob2 = ODEProblem(evolucion,rho0,tspan) 

sol2 = solve(prob2; save_everystep = true, saveat=0.01) 

peet = abs.(sol2[4,:])

plot!(sol2.t,peet, label ="Pe, Γ = 0.5π, δ = 0π")


omega = pi
gamma = pi
delta = 0*pi
ee = 0.0 + im*0.0 
ge = 0.0 + im*0.0 
eg = 0.0 + im*0.0 
gg = 1.0 + im*0.0 
rho0 = [gg,eg,ge,ee]
tspan =(0.0,4.0) 


function evolucion(drho,rho,p,t)
 drho[1] = -0.5*im*omega*(rho[2] - rho[3]) + gamma*rho[4] 
 drho[2] = 0.5*im*omega*(rho[4] - rho[1]) - 0.5*gamma*rho[2] + im*delta*rho[2]
 drho[3] = -0.5*im*omega*(rho[4] - rho[1]) - 0.5*gamma*rho[3] - im*delta*rho[3]
 drho[4] = 0.5*im*omega*(rho[2] - rho[3]) - gamma*rho[4] 
end 

prob2 = ODEProblem(evolucion,rho0,tspan) 

sol2 = solve(prob2; save_everystep = true, saveat=0.01) 

peet = abs.(sol2[4,:])

plot!(sol2.t,peet, label ="Pe, Γ = π, δ = 0π")


omega = pi
gamma = 2*pi
delta = 0*pi
ee = 0.0 + im*0.0 
ge = 0.0 + im*0.0 
eg = 0.0 + im*0.0 
gg = 1.0 + im*0.0 
rho0 = [gg,eg,ge,ee]
tspan =(0.0,4.0) 


function evolucion(drho,rho,p,t)
 drho[1] = -0.5*im*omega*(rho[2] - rho[3]) + gamma*rho[4] 
 drho[2] = 0.5*im*omega*(rho[4] - rho[1]) - 0.5*gamma*rho[2] + im*delta*rho[2]
 drho[3] = -0.5*im*omega*(rho[4] - rho[1]) - 0.5*gamma*rho[3] - im*delta*rho[3]
 drho[4] = 0.5*im*omega*(rho[2] - rho[3]) - gamma*rho[4] 
end 

prob2 = ODEProblem(evolucion,rho0,tspan) 

sol2 = solve(prob2; save_everystep = true, saveat=0.01) 

peet = abs.(sol2[4,:])
pggt = abs.(.*(sol2[1,:],conj(sol2[1,:])))
pegt = abs.(.*(sol2[2,:],conj(sol2[2,:])))
pget = abs.(.*(sol2[3,:],conj(sol2[3,:])))

plot!(sol2.t,peet, label ="Pe, Γ = 2π, δ = 0π")


omega = pi
gamma = 0.2*pi
delta = 3*pi
ee = 0.0 + im*0.0 
ge = 0.0 + im*0.0 
eg = 0.0 + im*0.0 
gg = 1.0 + im*0.0 
rho0 = [gg,eg,ge,ee]
tspan =(0.0,4.0) 


function evolucion(drho,rho,p,t)
 drho[1] = -0.5*im*omega*(rho[2] - rho[3]) + gamma*rho[4] 
 drho[2] = 0.5*im*omega*(rho[4] - rho[1]) - 0.5*gamma*rho[2] + im*delta*rho[2]
 drho[3] = -0.5*im*omega*(rho[4] - rho[1]) - 0.5*gamma*rho[3] - im*delta*rho[3]
 drho[4] = 0.5*im*omega*(rho[2] - rho[3]) - gamma*rho[4] 
end 

prob2 = ODEProblem(evolucion,rho0,tspan) 

sol2 = solve(prob2; save_everystep = true, saveat=0.01) 

peet = abs.(sol2[4,:])
pggt = abs.(.*(sol2[1,:],conj(sol2[1,:])))
pegt = abs.(.*(sol2[2,:],conj(sol2[2,:])))
pget = abs.(.*(sol2[3,:],conj(sol2[3,:])))

plot!(sol2.t,peet, label ="Pe, Γ = 0.2π, δ = 3π")











