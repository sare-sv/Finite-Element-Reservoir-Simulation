"""
Finite Element Simulation of Steady-State Single-Phase Flow
in a Heterogeneous 1D Oil Reservoir

Description:
This project implements linear and quadratic finite element methods
to simulate pressure distribution in a heterogeneous oil reservoir.
The model includes production wells, water injection, pressure
boundary conditions, and injection optimization based on bubble point pressure.

Author:
Sare Vahedi

Course:
Finite Element Method - Petroleum Engineering
"""

# Libraries:
import numpy as np
import matplotlib.pyplot as plt

# Input data:
# Reservoir geometry:
L = 9000.0
width = 1000.0
height = 100.0
A = width * height          # Cross-sectional flow area
# Fluid and rock properties:
Bo = 1.0
co = 3.5e-6
phi = 0.18
mu = 5.0
# Reservoir permeability:
k_zone1 = 60.0
k_zone2 = 20.0
x_interface = 4500.0
# Well locations:
x_well1 = 1500.0            # Producer 1
x_well_inj = 4500.0         # Injector
x_well2 = 7500.0            # Producer 2
# Well operating conditions:
q_prod = 1500
# Boundary conditions:
P_left = 4000.0
P_right = 4000.0
P_bubble = 2200.0
# Field-unit conversion factor:
darcy_const = 1.127e-3


# Mesh size:
# Number of linear elements:
n_elem = int(input("Enter the number of linear elements: "))

print("Number of elements used:", n_elem)


# Linear mesh:
n_node_lin = n_elem + 1
h_lin = L / n_elem
node_x_lin = np.zeros(n_node_lin)

for i in range(n_node_lin):
    node_x_lin[i] = i * h_lin

# Well node indices:
i_well1_lin = round(x_well1 / h_lin)
i_well_inj_lin = round(x_well_inj / h_lin)
i_well2_lin = round(x_well2 / h_lin)


# Linear stiffness matrix assembly:
K_lin = np.zeros((n_node_lin, n_node_lin))

for e in range(n_elem):
    x_left = node_x_lin[e]
    x_right = node_x_lin[e + 1]
    x_mid = (x_left + x_right) / 2.0

    # Element permeability:
    if x_mid < x_interface:
        k_e = k_zone1
    else:
        k_e = k_zone2

    C_e = darcy_const * k_e * A / (mu * Bo)

    # Local stiffness matrix:
    k_local = (C_e / h_lin) * np.array([[1.0, -1.0],
                                        [-1.0, 1.0]])

    # Assemble the global matrix
    K_lin[e, e] = K_lin[e, e] + k_local[0, 0]
    K_lin[e, e + 1] = K_lin[e, e + 1] + k_local[0, 1]
    K_lin[e + 1, e] = K_lin[e + 1, e] + k_local[1, 0]
    K_lin[e + 1, e + 1] = K_lin[e + 1, e + 1] + k_local[1, 1]
K_old = K_lin.copy()


# Initial solution but without injection
# Load vector:
F_no_inj = np.zeros(n_node_lin)
F_no_inj[i_well1_lin] = F_no_inj[i_well1_lin] - q_prod
F_no_inj[i_well2_lin] = F_no_inj[i_well2_lin] - q_prod
F_temp = F_no_inj.copy()
K_temp = K_old.copy()

# Apply boundary conditions:
F_temp = F_temp - K_temp[:, 0] * P_left
K_temp[:, 0] = 0.0
K_temp[0, :] = 0.0
K_temp[0, 0] = 1.0
F_temp[0] = P_left
F_temp = F_temp - K_temp[:, n_node_lin - 1] * P_right
K_temp[:, n_node_lin - 1] = 0.0
K_temp[n_node_lin - 1, :] = 0.0
K_temp[n_node_lin - 1, n_node_lin - 1] = 1.0
F_temp[n_node_lin - 1] = P_right

P_no_inj = np.linalg.solve(K_temp, F_temp)

# Check if pressure falls below the bubble point: 
min_p_no_inj = P_no_inj.min()

print("\nInitial solution without water injection")
print("Minimum reservoir pressure:", round(min_p_no_inj, 2), "psi")

if min_p_no_inj < P_bubble:
    print("Pressure is below the bubble point.")
    print("Water injection is required.")
else:
    print("Pressure is above the bubble point.")
    print("Water injection is not required.")


# Determine the minimum injection rate by increasing the rate, until the pressure is above the bubble point.

q_inj_try = 0.0
step = 1.0
while True:

    # Load vector
    F_try = np.zeros(n_node_lin)
    F_try[i_well1_lin] -= q_prod
    F_try[i_well2_lin] -= q_prod
    F_try[i_well_inj_lin] += q_inj_try

    F_temp = F_try.copy()
    K_temp = K_old.copy()

    # Apply boundary conditions
    F_temp = F_temp - K_temp[:, 0] * P_left
    K_temp[:, 0] = 0.0
    K_temp[0, :] = 0.0
    K_temp[0, 0] = 1.0
    F_temp[0] = P_left

    F_temp = F_temp - K_temp[:, n_node_lin - 1] * P_right
    K_temp[:, n_node_lin - 1] = 0.0
    K_temp[n_node_lin - 1, :] = 0.0
    K_temp[n_node_lin - 1, n_node_lin - 1] = 1.0
    F_temp[n_node_lin - 1] = P_right

    # Solve for pressure:
    P_try = np.linalg.solve(K_temp, F_temp)

    # Stop when the minimum pressure reaches the bubble point:
    if P_try.min() >= P_bubble:
        break

    q_inj_try += step

    # Safety check:
    if q_inj_try > 20000:
        print("Error: Target pressure was not reached.")
        break
q_inj_min = q_inj_try

print("\nMinimum water injection rate:")
print("q_inj_min =", round(q_inj_min, 1), "STB/day")


# Final linear FEM solution with optimal injection:
# Load vector
F_lin = np.zeros(n_node_lin)
F_lin[i_well1_lin] -= q_prod
F_lin[i_well2_lin] -= q_prod
F_lin[i_well_inj_lin] += q_inj_min
F_temp = F_lin.copy()
K_temp = K_old.copy()

# Apply boundary conditions:
F_temp = F_temp - K_temp[:, 0] * P_left
K_temp[:, 0] = 0.0
K_temp[0, :] = 0.0
K_temp[0, 0] = 1.0
F_temp[0] = P_left
F_temp = F_temp - K_temp[:, n_node_lin - 1] * P_right
K_temp[:, n_node_lin - 1] = 0.0
K_temp[n_node_lin - 1, :] = 0.0
K_temp[n_node_lin - 1, n_node_lin - 1] = 1.0
F_temp[n_node_lin - 1] = P_right

P_lin = np.linalg.solve(K_temp, F_temp)

# Calculate boundary water influx
R_lin = K_old.dot(P_lin) - F_lin
influx_left_lin = R_lin[0]
influx_right_lin = R_lin[-1]

print("\nFinal Linear FEM solution")
print("Minimum pressure:", round(P_lin.min(), 2), "psi")

print("\nBoundary water influx:")
print("Left boundary :", round(influx_left_lin, 2), "STB/day")
print("Right boundary:", round(influx_right_lin, 2), "STB/day")

# MB check:
total_in_lin = influx_left_lin + influx_right_lin + q_inj_min
total_out_lin = 2.0 * q_prod

print("\nMaterial balance check:")
print("Total inflow :", round(total_in_lin, 2), "STB/day")
print("Total outflow:", round(total_out_lin, 2), "STB/day")


# quadratic FEM (3 node elements)
# Each element has three nodes: left, middle and right, so the total number of nodes is 2*n_elem + 1.
n_node_quad = 2 * n_elem + 1
h_quad_elem = L / n_elem
node_x_quad = np.zeros(n_node_quad)

for i in range(n_node_quad):
    node_x_quad[i] = i * (h_quad_elem / 2.0)

i_well1_quad = round(x_well1 / (h_quad_elem / 2.0))
i_well_inj_quad = round(x_well_inj / (h_quad_elem / 2.0))
i_well2_quad = round(x_well2 / (h_quad_elem / 2.0))
K_quad = np.zeros((n_node_quad, n_node_quad))

for e in range(n_elem):

    node1 = 2 * e
    node2 = 2 * e + 1
    node3 = 2 * e + 2
    x_left = node_x_quad[node1]
    x_right = node_x_quad[node3]
    x_mid = (x_left + x_right) / 2.0

    # Determine element permeability:
    if x_mid < x_interface:
        k_e = k_zone1
    else:
        k_e = k_zone2
    C_e = darcy_const * k_e * A / (mu * Bo)

    # Local stiffness matrix:
    k_local_q = (C_e / (3.0 * h_quad_elem)) * np.array([
        [7.0, -8.0, 1.0],
        [-8.0, 16.0, -8.0],
        [1.0, -8.0, 7.0]
    ])
    idx = [node1, node2, node3]

    for a in range(3):
        for b in range(3):
            K_quad[idx[a], idx[b]] += k_local_q[a, b]
K_quad_old = K_quad.copy()

# Solve the quadratic FEM system:
F_quad = np.zeros(n_node_quad)
F_quad[i_well1_quad] -= q_prod
F_quad[i_well2_quad] -= q_prod
F_quad[i_well_inj_quad] += q_inj_min
F_temp_q = F_quad.copy()
K_temp_q = K_quad_old.copy()

# Apply boundary conditions:
F_temp_q = F_temp_q - K_temp_q[:, 0] * P_left
K_temp_q[:, 0] = 0.0
K_temp_q[0, :] = 0.0
K_temp_q[0, 0] = 1.0
F_temp_q[0] = P_left
last = n_node_quad - 1
F_temp_q = F_temp_q - K_temp_q[:, last] * P_right
K_temp_q[:, last] = 0.0
K_temp_q[last, :] = 0.0
K_temp_q[last, last] = 1.0
F_temp_q[last] = P_right

P_quad = np.linalg.solve(K_temp_q, F_temp_q)

# Calculate boundary water influx:
R_quad = K_quad_old.dot(P_quad) - F_quad
influx_left_quad = R_quad[0]
influx_right_quad = R_quad[-1]

print("\nQuadratic FEM results:")
print("Left boundary influx :", round(influx_left_quad, 2), "STB/day")
print("Right boundary influx:", round(influx_right_quad, 2), "STB/day")


# Compare Linear and Quadratic FEM Results:
print("\nPressure comparison (Linear vs. Quadratic):")
print("   x (ft)      Linear (psi)      Quadratic (psi)      Difference")

for i in range(n_node_lin):
    j = 2 * i
    diff = P_lin[i] - P_quad[j]
    print("  ", round(node_x_lin[i], 1),
          "     ", round(P_lin[i], 3),
          "     ", round(P_quad[j], 3),
          "     ", round(diff, 5))

print("\nThe two solutions are almost identical at the shared nodes.")
print("This is expected because the pressure distribution is nearly linear in each zone.")

# Pressure distribution Plot:
fig, ax = plt.subplots(figsize=(9, 5.5))

# Plot pressure distribution:
ax.plot(node_x_lin, P_lin, "o-", color="tab:blue",
        label="Linear Element", linewidth=2, markersize=6)
ax.plot(node_x_quad, P_quad, "s--", color="tab:orange",
        label="Quadratic Element", linewidth=1.5, markersize=4)

# Bubble point pressure:
ax.axhline(y=P_bubble, color="red", linestyle=":",
           linewidth=2, label="Bubble Point (2200 psi)")

# Well locations:
ax.axvline(x=x_well1, color="green", linestyle="--", alpha=0.5)
ax.axvline(x=x_well_inj, color="purple", linestyle="--", alpha=0.5)
ax.axvline(x=x_well2, color="green", linestyle="--", alpha=0.5)

y_bottom, y_top = ax.get_ylim()
ax.set_ylim(y_bottom, y_top + 0.12 * (y_top - y_bottom)) 
ax.text(x_well1, y_top, "Producer 1", ha="center", va="bottom", fontsize=9)
ax.text(x_well_inj, y_top, "Injector", ha="center", va="bottom", fontsize=9)
ax.text(x_well2, y_top, "Producer 2", ha="center", va="bottom", fontsize=9)

# Plot settings:
ax.set_xlabel("Reservoir Length (ft)")
ax.set_ylabel("Pressure (psi)")
ax.set_title("Pressure Distribution")
ax.legend()
ax.grid(True, alpha=0.3)
plt.tight_layout()

# Save and display the figure
plt.savefig("reservoir_pressure_plot.png", dpi=150)
plt.show()
print("\nPressure plot saved as 'reservoir_pressure_plot.png'.")

# Results Summary
print("")
print("FINAL RESULTS SUMMARY")
print("")
print("Number of elements used:", n_elem)
print("Minimum water injection rate required:", round(q_inj_min, 1), "STB/day")
print("Minimum reservoir pressure after injection:", round(P_lin.min(), 2), "psi")
print("Water influx at left boundary:", round(influx_left_lin, 2), "STB/day")
print("Water influx at right boundary:", round(influx_right_lin, 2), "STB/day")
