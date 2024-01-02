using Plots
using PlotUtils
pythonplot()

r_min = 1
r_max = 10
m = 100
k = 300
dt = 0.0001
playback_speed = 1

n = Int(ceil(log(m/(m - 2*pi), r_max/r_min)))
r_max = r_min*(m/(m - 2*pi))^n

c_v = 717
R = 287

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
x = zeros(m, n)
y = zeros(m, n)

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

                drhodxi = rho[i+1, j, t] - rho[i, j, t]
                dudxi = u[i+1, j, t] - u[i, j, t]
                dvdxi = v[i+1, j, t] - v[i, j, t]
                dedxi = e[i+1, j, t] - e[i, j, t]
                dpdxi = p[i+1, j, t] - p[i, j, t]

            elseif i == 1

                drhodxi = rho[i+1, j, t] - rho[i, j, t]
                dudxi = u[i+1, j, t] - u[i, j, t]
                dvdxi = v[i+1, j, t] - v[i, j, t]
                dedxi = e[i+1, j, t] - e[i, j, t]
                dpdxi = p[i+1, j, t] - p[i, j, t]

            elseif i == m

                drhodxi = rho[1, j, t] - rho[i, j, t]
                dudxi = u[1, j, t] - u[i, j, t]
                dvdxi = v[1, j, t] - v[i, j, t]
                dedxi = e[1, j, t] - e[i, j, t]
                dpdxi = p[1, j, t] - p[i, j, t]

            end

            if j != 1 && j != n

                drhodeta = rho[i, j+1, t] - rho[i, j, t]
                dudeta = u[i, j+1, t] - u[i, j, t]
                dvdeta = v[i, j+1, t] - v[i, j, t]
                dedeta = e[i, j+1, t] - e[i, j, t]
                dpdeta = p[i, j+1, t] - p[i, j, t]
                        
            elseif j == 1

                drhodeta = (-3*rho[i, j, t] + 4*rho[i, j+1, t] - rho[i, j+2, t])/2
                dudeta = (-3*u[i, j, t] + 4*u[i, j+1, t] - u[i, j+2, t])/2
                dvdeta = (-3*v[i, j, t] + 4*v[i, j+1, t] - v[i, j+2, t])/2
                dedeta = (-3*e[i, j, t] + 4*e[i, j+1, t] - e[i, j+2, t])/2
                dpdeta = 0

            elseif j == n

                drhodeta = (3*rho[i, j, t] - 4*rho[i, j-1, t] + rho[i, j-2, t])/2
                dudeta = (3*u[i, j, t] - 4*u[i, j-1, t] + u[i, j-2, t])/2
                dvdeta = (3*v[i, j, t] - 4*v[i, j-1, t] + v[i, j-2, t])/2
                dedeta = (3*e[i, j, t] - 4*e[i, j-1, t] + e[i, j-2, t])/2
                dpdeta = 0

            end

            xi = i - 1
            eta = j - 1

            dxdxi = -(2*r_min*pi*sin(2*pi*xi/m)*(m/(m - 2*pi))^eta)/m
            dydxi = (2*r_min*pi*cos(2*pi*xi/m)*(m/(m - 2*pi))^eta)/m
            dxdeta = r_min*log(m/(m - 2*pi))*cos(2*pi*xi/m)*(m/(m - 2*pi))^eta
            dydeta = r_min*log(m/(m - 2*pi))*sin(2*pi*xi/m)*(m/(m - 2*pi))^eta

            J = -(2*r_min^2*pi*log(m/(m - 2*pi))*(m/(m - 2*pi))^(2*eta))/m

            drhodx = 1/J*(drhodxi*dydeta - drhodeta*dydxi)
            drhody = 1/J*(drhodeta*dxdxi - drhodxi*dxdeta)
            dudx = 1/J*(dudxi*dydeta - dudeta*dydxi)
            dudy = 1/J*(dudeta*dxdxi - dudxi*dxdeta)
            dvdx = 1/J*(dvdxi*dydeta - dvdeta*dydxi)
            dvdy = 1/J*(dvdeta*dxdxi - dvdxi*dxdeta)
            dedx = 1/J*(dedxi*dydeta - dedeta*dydxi)
            dedy = 1/J*(dedeta*dxdxi - dedxi*dxdeta)
            dpdx = 1/J*(dpdxi*dydeta - dpdeta*dydxi)
            dpdy = 1/J*(dpdeta*dxdxi - dpdxi*dxdeta)

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

                drhodxi = rho_bar[i, j] - rho_bar[i-1, j]
                dudxi = u_bar[i, j] - u_bar[i-1, j]
                dvdxi = v_bar[i, j] - v_bar[i-1, j]
                dedxi = e_bar[i, j] - e_bar[i-1, j]
                dpdxi = p_bar[i, j] - p_bar[i-1, j]

            elseif i == 1

                drhodxi = rho_bar[i, j] - rho_bar[m, j]
                dudxi = u_bar[i, j] - u_bar[m, j]
                dvdxi = v_bar[i, j] - v_bar[m, j]
                dedxi = e_bar[i, j] - e_bar[m, j]
                dpdxi = p_bar[i, j] - p_bar[m, j]

            elseif i == m

                drhodxi = rho_bar[i, j] - rho_bar[i-1, j]
                dudxi = u_bar[i, j] - u_bar[i-1, j]
                dvdxi = v_bar[i, j] - v_bar[i-1, j]
                dedxi = e_bar[i, j] - e_bar[i-1, j]
                dpdxi = p_bar[i, j] - p_bar[i-1, j]

            end

            if j != 1 && j != n

                drhodeta = rho_bar[i, j] - rho_bar[i, j-1]
                dudeta = u_bar[i, j] - u_bar[i, j-1]
                dvdeta = v_bar[i, j] - v_bar[i, j-1]
                dedeta = e_bar[i, j] - e_bar[i, j-1]
                dpdeta = p_bar[i, j] - p_bar[i, j-1]
                        
            elseif j == 1

                drhodeta = (-3*rho_bar[i, j] + 4*rho_bar[i, j+1] - rho_bar[i, j+2])/2
                dudeta = (-3*u_bar[i, j] + 4*u_bar[i, j+1] - u_bar[i, j+2])/2
                dvdeta = (-3*v_bar[i, j] + 4*v_bar[i, j+1] - v_bar[i, j+2])/2
                dedeta = (-3*e_bar[i, j] + 4*e_bar[i, j+1] - e_bar[i, j+2])/2
                dpdeta = 0

            elseif j == n

                drhodeta = (3*rho_bar[i, j] - 4*rho_bar[i, j-1] + rho_bar[i, j-2])/2
                dudeta = (3*u_bar[i, j] - 4*u_bar[i, j-1] + u_bar[i, j-2])/2
                dvdeta = (3*v_bar[i, j] - 4*v_bar[i, j-1] + v_bar[i, j-2])/2
                dedeta = (3*e_bar[i, j] - 4*e_bar[i, j-1] + e_bar[i, j-2])/2
                dpdeta = 0

            end

            xi = i - 1
            eta = j - 1

            dxdxi = -(2*r_min*pi*sin(2*pi*xi/m)*(m/(m - 2*pi))^eta)/m
            dydxi = (2*r_min*pi*cos(2*pi*xi/m)*(m/(m - 2*pi))^eta)/m
            dxdeta = r_min*log(m/(m - 2*pi))*cos(2*pi*xi/m)*(m/(m - 2*pi))^eta
            dydeta = r_min*log(m/(m - 2*pi))*sin(2*pi*xi/m)*(m/(m - 2*pi))^eta

            J = -(2*r_min^2*pi*log(m/(m - 2*pi))*(m/(m - 2*pi))^(2*eta))/m

            drhodx = 1/J*(drhodxi*dydeta - drhodeta*dydxi)
            drhody = 1/J*(drhodeta*dxdxi - drhodxi*dxdeta)
            dudx = 1/J*(dudxi*dydeta - dudeta*dydxi)
            dudy = 1/J*(dudeta*dxdxi - dudxi*dxdeta)
            dvdx = 1/J*(dvdxi*dydeta - dvdeta*dydxi)
            dvdy = 1/J*(dvdeta*dxdxi - dvdxi*dxdeta)
            dedx = 1/J*(dedxi*dydeta - dedeta*dydxi)
            dedy = 1/J*(dedeta*dxdxi - dedxi*dxdeta)
            dpdx = 1/J*(dpdxi*dydeta - dpdeta*dydxi)
            dpdy = 1/J*(dpdeta*dxdxi - dpdxi*dxdeta)

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

            if isnan(p[i, j, t+1]) == true
                error("NaN @ t = ", t)
            end

        end
    end

end

r = r_min*(m/(m - 2*pi)).^(0:n-1)
theta = range(0, 2*pi, m)

for i = eachindex(theta)
    for j = eachindex(r)

        x[i, j] = r[j]*cos(theta[i])
        y[i, j] = r[j]*sin(theta[i])

    end
end

p_min = Int(floor(minimum(p)))
p_max = Int(ceil(maximum(p)))

anim = @animate for i = 1:playback_speed:k
    contourf(x, y, p[:, :, i], color=:turbo, levels=range(p_min, p_max, 100), colorbar_ticks=(range(p_min, p_max, 10), round.(range(p_min, p_max, 10), digits=2)), aspect_ratio=:equal, clim=(p_min, p_max))
end

gif(anim, "CFD_transformation.gif")
