using Plots

m = 100
k = 300
dx = 0.1
dt = 0.0001

rho = zeros(m, k)
u = zeros(m, k)
e = zeros(m, k)
T = zeros(m, k)
p = zeros(m, k)
rho_bar = zeros(m)
u_bar = zeros(m)
e_bar = zeros(m)
T_bar = zeros(m)
p_bar = zeros(m)
drhodt = zeros(m)
dudt = zeros(m)
dedt = zeros(m)

c_v = 717
R = 287

x = range(-5, 5, m)

T[:, 1] .= 273.15
p[:, 1] = 10^5 .+ 10*exp.(-x.^2)
e[:, 1] = c_v*T[:, 1]
rho[:, 1] = p[:, 1]./(R*T[:, 1])
u[:, 1] .= 0

for t = 1:k-1

    for i = 1:m

        if i != 1 && i != m

            drhodx = (rho[i+1, t] - rho[i, t])/dx
            dudx = (u[i+1, t] - u[i, t])/dx
            dedx = (e[i+1, t] - e[i, t])/dx
            dpdx = (p[i+1, t] - p[i, t])/dx

        elseif i == 1

            drhodx = (-3*rho[i, t] + 4*rho[i+1, t] - rho[i+2, t])/(2*dx)
            dudx = (-3*u[i, t] + 4*u[i+1, t] - u[i+2, t])/(2*dx)
            dedx = (-3*e[i, t] + 4*e[i+1, t] - e[i+2, t])/(2*dx)
            dpdx = 0
            
        elseif i == m

            drhodx = (3*rho[i, t] - 4*rho[i-1, t] + rho[i-2, t])/(2*dx)
            dudx = (3*u[i, t] - 4*u[i-1, t] + u[i-2, t])/(2*dx)
            dedx = (3*e[i, t] - 4*e[i-1, t] + e[i-2, t])/(2*dx)
            dpdx = 0
        
        end

        drhodt[i] = -(rho[i, t]*dudx + u[i, t]*drhodx)
        dudt[i] = -(u[i, t]*dudx + 1/rho[i, t]*dpdx)
        dedt[i] = -(u[i, t]*dedx + p[i, t]/rho[i, t]*dudx)

        rho_bar[i] = rho[i, t] + drhodt[i]*dt
        u_bar[i] = u[i, t] + dudt[i]*dt
        e_bar[i] = e[i, t] + dedt[i]*dt
        T_bar[i] = e_bar[i]/c_v
        p_bar[i] = rho_bar[i]*R*T_bar[i]

    end

    for i = 1:m

        if i != 1 && i != m

            drhodx = (rho_bar[i] - rho_bar[i-1])/dx
            dudx = (u_bar[i] - u_bar[i-1])/dx
            dedx = (e_bar[i] - e_bar[i-1])/dx
            dpdx = (p_bar[i] - p_bar[i-1])/dx
        
        elseif i == 1

            drhodx = (-3*rho_bar[i] + 4*rho_bar[i+1] - rho_bar[i+2])/(2*dx)
            dudx = (-3*u_bar[i] + 4*u_bar[i+1] - u_bar[i+2])/(2*dx)
            dedx = (-3*e_bar[i] + 4*e_bar[i+1] - e_bar[i+2])/(2*dx)
            dpdx = 0
            
        elseif i == m

            drhodx = (3*rho_bar[i] - 4*rho_bar[i-1] + rho_bar[i-2])/(2*dx)
            dudx = (3*u_bar[i] - 4*u_bar[i-1] + u_bar[i-2])/(2*dx)
            dedx = (3*e_bar[i] - 4*e_bar[i-1] + e_bar[i-2])/(2*dx)
            dpdx = 0
        
        end

        drhodt_bar = -(rho_bar[i]*dudx + u_bar[i]*drhodx)
        dudt_bar = -(u_bar[i]*dudx + 1/rho_bar[i]*dpdx)
        dedt_bar = -(u_bar[i]*dedx + p_bar[i]/rho_bar[i]*dudx)

        drhodt_av = 1/2*(drhodt[i] + drhodt_bar)
        dudt_av = 1/2*(dudt[i] + dudt_bar)
        dedt_av = 1/2*(dedt[i] + dedt_bar)

        rho[i, t+1] = rho[i, t] + drhodt_av*dt
        u[i, t+1] = u[i, t] + dudt_av*dt
        e[i, t+1] = e[i, t] + dedt_av*dt
        T[i, t+1] = e[i, t+1]/c_v
        p[i, t+1] = rho[i, t+1]*R*T[i, t+1]

    end

end

anim = @animate for i = 1:k
    plot(p[:, i], ylim = (100000, 100010))
end

gif(anim, "CFD1D.gif", fps = 30)