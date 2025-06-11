import matplotlib.pyplot as plt
import numpy as np
import matplotlib.patches as pat

angulo = np.pi/4

def lineplot(x,y,z, thickness, style, marker):
    plt.plot((x - np.multiply(z, np.cos(angulo))), (y - np.multiply(z, np.sin(angulo))),color='black',linewidth=thickness,linestyle=style,marker=marker,)

def arrowplot(vector, length, thickness):
    plt.arrow((vector[0] - vector[2] * np.cos(angulo)), (vector[1] - vector[2] * np.sin(angulo)), (length[0] - length[2] * np.cos(angulo)), (length[1] - length[2] * np.sin(angulo)),color='black',width=thickness,head_width=6*thickness,length_includes_head=True)



# build image

fig = plt.figure(figsize=(10, 6))
ax = fig.add_subplot()

endx = 3.0
endy = 1.5
endz = 1.0

# plots axis
origin = [0, 0, 0]
arrowplot(origin, [3.5,0,0], 0.01)
arrowplot(origin, [0,2,0], 0.01)
arrowplot(origin, [0,0,1.5], 0.01)
origin_b = [endx,endy,endz]
arrowplot(origin_b, [1,0,0], 0.01)
arrowplot(origin_b, [0,1,0], 0.01)
arrowplot(origin_b, [0,0,1], 0.01)

# plots lines
lineplot([0,endx],[0,endy],[0,endz],1,'--','.')
lineplot([0,endx],[endy,endy],[0,endz],1,'--','.')
lineplot([endx,endx],[0,endy],[endz,endz],1,'--','.')
lineplot([endx,endx],[0,0],[0,endz],1,'--','.')
lineplot([0,endx],[0,0],[endz,endz],1,'--','.')
lineplot([0,endx],[0,0],[0,endz],1,'--','.')

# plots arcs
a1 = pat.Arc([0,0], 0.35, 0.38, angle=225, theta1=0, theta2=118)
a2 = pat.Arc([0,0], 0.35, 0.35, angle=20, theta1=0, theta2=70)
ax.add_patch(a1)
ax.add_patch(a2)

# plots texts
size = 'x-large'
plt.text(-0.22, 0, r"$\mathcal{O}_P$", fontsize=size)
plt.text(endx + 0.5, 0.1, "$y$", fontsize=size)
plt.text(endx - 0.05, 0.1, "$y_0$", fontsize=size)
plt.text(-0.15, endy + 0.5, "$z$", fontsize=size)
plt.text(-0.2, endy, "$z_0$", fontsize=size)
plt.text(-1.2, -1.05, "$x$", fontsize=size)
plt.text(-0.95, -0.75, "$x_0$", fontsize=size)
plt.text(0.05, -0.3, "$\phi_0$", fontsize=size)
plt.text(0.1, 0.2, "$\\theta_0$", fontsize=size)
plt.text(endx - 0.65, endy - 0.9, r"$\mathcal{O}_B$", fontsize=size)
plt.text(0.9, -0.45, "$\\rho_0$", fontsize=size)
plt.text(endx - 1.55, endy - 1.4, "$u$", fontsize=size)
plt.text(endx + 0.25, endy - 0.65, "$v$", fontsize=size)
plt.text(endx - 0.85, endy + 0.3, "$w$", fontsize=size)
plt.text(endx - 1.3, endy - 2.4, "$(x_0,y_0) \\equiv (\\rho_0,\phi_0)$", fontsize=size)
plt.text(0.9, 0.45, "$r_0$", fontsize=size)

plt.axis('equal')
plt.axis('off')
plt.savefig('sistema_coordenadas', dpi=800)