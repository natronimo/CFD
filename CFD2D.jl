using Plots

m = 100
n = 100
k = 300
dx = 0.1
dy = 0.1
dt = 0.0001

rho = zeros(m, n, k)
u = zeros(m, n, k)
v = zeros(m, n, k)
e = zeros(m, n, k)
T = zeros(m, n, k)
p = zeros(m, n, k)
rho_bar = zeros(m, n)
u_bar = zeros(m, n)
v_bar = zeros(m, n)
e_bar = zeros(m, n)
T_bar = zeros(m, n)
p_bar = zeros(m, n)
drhodt = zeros(m, n)
dudt = zeros(m, n)
dvdt = zeros(m, n)
dedt = zeros(m, n)

c_v = 717
R = 287

T[:, :, 1] .= 273.15
for i = 1:m
    x = 10*i/m - 5
    for j = 1:n
        y = 10*j/n - 5
        p[i, j, 1] = 10^5 + 10*exp(-(x^2 + y^2))
    end
end
e[:, :, 1] = c_v*T[:, :, 1]
rho[:, :, 1] = p[:, :, 1]./(R*T[:, :, 1])
u[:, :, 1] .= 0

for t = 1:k-1

    for i = 1:m
        for j = 1:n

            if i != 1 && i != m

                drhodx = (rho[i+1, j, t] - rho[i, j, t])/dx
                dudx = (u[i+1, j, t] - u[i, j, t])/dx
                dvdx = (v[i+1, j, t] - v[i, j, t])/dx
                dedx = (e[i+1, j, t] - e[i, j, t])/dx
                dpdx = (p[i+1, j, t] - p[i, j, t])/dx

            elseif i == 1

                drhodx = (-3*rho[i, j, t] + 4*rho[i+1, j, t] - rho[i+2, j, t])/(2*dx)
                dudx = (-3*u[i, j, t] + 4*u[i+1, j, t] - u[i+2, j, t])/(2*dx)
                dvdx = (-3*v[i, j, t] + 4*v[i+1, j, t] - v[i+2, j, t])/(2*dx)
                dedx = (-3*e[i, j, t] + 4*e[i+1, j, t] - e[i+2, j, t])/(2*dx)
                dpdx = 0

            elseif i == m

                drhodx = (3*rho[i, j, t] - 4*rho[i-1, j, t] + rho[i-2, j, t])/(2*dx)
                dudx = (3*u[i, j, t] - 4*u[i-1, j, t] + u[i-2, j, t])/(2*dx)
                dvdx = (3*v[i, j, t] - 4*v[i-1, j, t] + v[i-2, j, t])/(2*dx)
                dedx = (3*e[i, j, t] - 4*e[i-1, j, t] + e[i-2, j, t])/(2*dx)
                dpdx = 0

            end

            if j != 1 && j != n

                drhody = (rho[i, j+1, t] - rho[i, j, t])/dy
                dudy = (u[i, j+1, t] - u[i, j, t])/dy
                dvdy = (v[i, j+1, t] - v[i, j, t])/dy
                dedy = (e[i, j+1, t] - e[i, j, t])/dy
                dpdy = (p[i, j+1, t] - p[i, j, t])/dy
                        
            elseif j == 1

                drhody = (-3*rho[i, j, t] + 4*rho[i, j+1, t] - rho[i, j+2, t])/(2*dy)
                dudy = (-3*u[i, j, t] + 4*u[i, j+1, t] - u[i, j+2, t])/(2*dy)
                dvdy = (-3*v[i, j, t] + 4*v[i, j+1, t] - v[i, j+2, t])/(2*dy)
                dedy = (-3*e[i, j, t] + 4*e[i, j+1, t] - e[i, j+2, t])/(2*dy)
                dpdy = 0

            elseif j == n

                drhody = (3*rho[i, j, t] - 4*rho[i, j-1, t] + rho[i, j-2, t])/(2*dy)
                dudy = (3*u[i, j, t] - 4*u[i, j-1, t] + u[i, j-2, t])/(2*dy)
                dvdy = (3*v[i, j, t] - 4*v[i, j-1, t] + v[i, j-2, t])/(2*dy)
                dedy = (3*e[i, j, t] - 4*e[i, j-1, t] + e[i, j-2, t])/(2*dy)
                dpdy = 0

            end

            drhodt[i, j] = -(rho[i, j, t]*dudx + u[i, j, t]*drhodx + rho[i, j, t]*dvdy + v[i, j, t]*drhody)
            dudt[i, j] = -(u[i, j, t]*dudx + v[i, j, t]*dudy + 1/rho[i, j, t]*dpdx)
            dvdt[i, j] = -(u[i, j, t]*dvdx + v[i, j, t]*dvdy + 1/rho[i, j, t]*dpdy)
            dedt[i, j] = -(u[i, j, t]*dedx + v[i, j, t]*dedy + p[i, j, t]/rho[i, j, t]*dudx + p[i, j, t]/rho[i, j, t]*dvdy)

            rho_bar[i, j] = rho[i, j, t] + drhodt[i, j]*dt
            u_bar[i, j] = u[i, j, t] + dudt[i, j]*dt
            v_bar[i, j] = v[i, j, t] + dvdt[i, j]*dt
            e_bar[i, j] = e[i, j, t] + dedt[i, j]*dt
            T_bar[i, j] = e_bar[i, j]/c_v
            p_bar[i, j] = rho_bar[i, j]*R*T_bar[i, j]

        end
    end

    for i = 1:m
        for j = 1:n

            if i != 1 && i != m

                drhodx = (rho_bar[i, j] - rho_bar[i-1, j])/dx
                dudx = (u_bar[i, j] - u_bar[i-1, j])/dx
                dvdx = (v_bar[i, j] - v_bar[i-1, j])/dx
                dedx = (e_bar[i, j] - e_bar[i-1, j])/dx
                dpdx = (p_bar[i, j] - p_bar[i-1, j])/dx

            elseif i == 1

                drhodx = (-3*rho_bar[i, j] + 4*rho_bar[i+1, j] - rho_bar[i+2, j])/(2*dx)
                dudx = (-3*u_bar[i, j] + 4*u_bar[i+1, j] - u_bar[i+2, j])/(2*dx)
                dvdx = (-3*v_bar[i, j] + 4*v_bar[i+1, j] - v_bar[i+2, j])/(2*dx)
                dedx = (-3*e_bar[i, j] + 4*e_bar[i+1, j] - e_bar[i+2, j])/(2*dx)
                dpdx = 0

            elseif i == m

                drhodx = (3*rho_bar[i, j] - 4*rho_bar[i-1, j] + rho_bar[i-2, j])/(2*dx)
                dudx = (3*u_bar[i, j] - 4*u_bar[i-1, j] + u_bar[i-2, j])/(2*dx)
                dvdx = (3*v_bar[i, j] - 4*v_bar[i-1, j] + v_bar[i-2, j])/(2*dx)
                dedx = (3*e_bar[i, j] - 4*e_bar[i-1, j] + e_bar[i-2, j])/(2*dx)
                dpdx = 0

            end

            if j != 1 && j != n

                drhody = (rho_bar[i, j] - rho_bar[i, j-1])/dy
                dudy = (u_bar[i, j] - u_bar[i, j-1])/dy
                dvdy = (v_bar[i, j] - v_bar[i, j-1])/dy
                dedy = (e_bar[i, j] - e_bar[i, j-1])/dy
                dpdy = (p_bar[i, j] - p_bar[i, j-1])/dy

            elseif j == 1

                drhody = (-3*rho_bar[i, j] + 4*rho_bar[i, j+1] - rho_bar[i, j+2])/(2*dy)
                dudy = (-3*u_bar[i, j] + 4*u_bar[i, j+1] - u_bar[i, j+2])/(2*dy)
                dvdy = (-3*v_bar[i, j] + 4*v_bar[i, j+1] - v_bar[i, j+2])/(2*dy)
                dedy = (-3*e_bar[i, j] + 4*e_bar[i, j+1] - e_bar[i, j+2])/(2*dy)
                dpdy = 0

            elseif j == n

                drhody = (3*rho_bar[i, j] - 4*rho_bar[i, j-1] + rho_bar[i, j-2])/(2*dy)
                dudy = (3*u_bar[i, j] - 4*u_bar[i, j-1] + u_bar[i, j-2])/(2*dy)
                dvdy = (3*v_bar[i, j] - 4*v_bar[i, j-1] + v_bar[i, j-2])/(2*dy)
                dedy = (3*e_bar[i, j] - 4*e_bar[i, j-1] + e_bar[i, j-2])/(2*dy)
                dpdy = 0

            end

            drhodt_bar = -(rho_bar[i, j]*dudx + u_bar[i, j]*drhodx + rho_bar[i, j]*dvdy + v_bar[i, j]*drhody)
            dudt_bar = -(u_bar[i, j]*dudx + v_bar[i, j]*dudy + 1/rho_bar[i, j]*dpdx)
            dvdt_bar = -(u_bar[i, j]*dvdx + v_bar[i, j]*dvdy + 1/rho_bar[i, j]*dpdy)
            dedt_bar = -(u_bar[i, j]*dedx + v_bar[i, j]*dedy + p_bar[i, j]/rho_bar[i, j]*dudx + p_bar[i, j]/rho_bar[i, j]*dvdy)

            drhodt_av = 1/2*(drhodt[i, j] + drhodt_bar)
            dudt_av = 1/2*(dudt[i, j] + dudt_bar)
            dvdt_av = 1/2*(dvdt[i, j] + dvdt_bar)
            dedt_av = 1/2*(dedt[i, j] + dedt_bar)

            rho[i, j, t+1] = rho[i, j, t] + drhodt_av*dt
            u[i, j, t+1] = u[i, j, t] + dudt_av*dt
            v[i, j, t+1] = v[i, j, t] + dvdt_av*dt
            e[i, j, t+1] = e[i, j, t] + dedt_av*dt
            T[i, j, t+1] = e[i, j, t+1]/c_v
            p[i, j, t+1] = rho[i, j, t+1]*R*T[i, j, t+1]

        end
    end

end

anim = @animate for i = 1:k
    surface(p[:, :, i], zlim = (minimum(p), maximum(p)))
end

gif(anim, "CFD2D.gif", fps = 30)