#!/usr/bin/env python3
import numpy as np
import sympy as sp



class forward_kinematic:

    def __init__(self, method):

        self.method = method

    def dh_transform_original(self, theta, d, a, alpha):
        """
        Computes the homogeneous transformation matrix using DH parameters.
        
        Parameters:
        theta (float): Joint angle (radians)
        d (float): Link offset
        a (float): Link length
        alpha (float): Link twist (radians)
        
        Returns:
        sympy.Matrix: 4x4 homogeneous transformation matrix
        """
        # Create the transformation matrix
        T = sp.Matrix([
            [sp.cos(theta), -sp.sin(theta)*sp.cos(alpha),  sp.sin(theta)*sp.sin(alpha), a*sp.cos(theta)],
            [sp.sin(theta),  sp.cos(theta)*sp.cos(alpha), -sp.cos(theta)*sp.sin(alpha), a*sp.sin(theta)],
            [0,              sp.sin(alpha),               sp.cos(alpha),               d],
            [0,              0,                           0,                           1]
        ])
        
        return sp.simplify(T)

    def dh_transformation_modified(self, theta, d, a, alpha):
        """
        Computes the homogeneous transformation matrix using DH parameters.
        
        Parameters:
        theta (float): Joint angle (radians)
        d (float): Link offset
        a (float): Link length
        alpha (float): Link twist (radians)
        
        Returns:
        sympy.Matrix: 4x4 homogeneous transformation matrix
        """
        # Create the transformation matrix
        T = sp.Matrix([
            [sp.cos(theta), -sp.sin(theta),  0, a],
            [sp.sin(theta)*sp.cos(alpha), sp.cos(theta)*sp.cos(alpha), -sp.sin(alpha), -d*sp.sin(alpha)],
            [sp.sin(theta)*sp.sin(alpha), sp.cos(theta)*sp.sin(alpha),  sp.cos(alpha),  d*sp.cos(alpha)],
            [0, 0, 0, 1]
        ])
        
        return sp.simplify(T)
    
    def skew_symmatric(self, v):
        """
        Returns the skew-symmetric matrix of a vector.
        
        Parameters:
        v (sympy.Matrix): 3x1 vector
        
        Returns:
        sympy.Matrix: 3x3 skew-symmetric matrix
        """
        return sp.Matrix([
            [0, -v[2], v[1]],
            [v[2], 0, -v[0]],
            [-v[1], v[0], 0]
        ])

    def twist_to_transform(self, twist, theta):
        """
        Returns the exponential of a twist vector.
        
        Parameters:
        twist (sympy.Matrix): 6x1 twist vector
        theta (float): Joint angle (radians)
        
        Returns:
        sympy.Matrix: 4x4 homogeneous transformation matrix
        """
        # Extract the angular and linear velocity components
        omega = twist[:3]
        v = twist[3:]
        
        # Compute the skew-symmetric matrix of the angular velocity
        omega_hat = self.skew_symmatric(omega)
        
        # Compute the rotation matrix
        R = sp.eye(3) + sp.sin(theta)*omega_hat + (1 - sp.cos(theta))*omega_hat**2
        
        # Compute the translation matrix
        p = (sp.eye(3)*theta + (1 - sp.cos(theta))*omega_hat + (theta - sp.sin(theta))*omega_hat**2)*v
        
        # Combine the rotation and translation matrices
        T = sp.eye(4)
        T[:3, :3] = R
        T[:3, 3] = p
        
        return sp.simplify(T)



    def transform_frame_to_end_effector(self, dh_table):
        
        if self.method == 1:

            T = sp.eye(4)  # Initialize as identity matrix
            for i in range(dh_table.shape[0]):
                theta, d, a, alpha = dh_table[i]
                T_i = self.dh_transform_original(theta, d, a, alpha)
                T = sp.simplify(T * T_i)  # Simplify at each step to reduce complexity
            
            return T
        
        elif self.method == 2:

            T = sp.eye(4)
            for i in range(dh_table.shape[0]):
                theta, d, a, alpha = dh_table[i]
                T_i = self.dh_transformation_modified(theta, d, a, alpha)
                T = sp.simplify(T * T_i)

        elif self.method == 3:
            



test_me = forward_kinematic(1)

# Define symbolic variables for DH parameters
l1, l2, l3, l4, l5, l6, l7 = sp.symbols('l1 l2 l3 l4 l5 l6 l7')
zeta_1, zeta_2, zeta_3, zeta_4 = sp.symbols('zeta_1 zeta_2 zeta_3 zeta_4')

# Define the DH table as a numpy array
dh_table = np.array([
    [0, sp.pi/2, 0, sp.pi/2],
    [l1, 0, 0, zeta_1],
    [0, -sp.pi/2, -l2, -sp.pi/2],
    [l3, 0, 0, zeta_2 + sp.pi/2],
    [l4, 0, l5, -sp.pi/2],
    [0, sp.pi/2, 0, zeta_3],
    [l6, 0, 0, zeta_4],
    [l7, 0, 0, 0]
])

# Compute the transformation matrix from base frame to end effector
T_end_effector = test_me.transform_frame_to_end_effector(dh_table)

# Simplify the x, y, z positions
x_position = sp.simplify(T_end_effector[0, 3])
y_position = sp.simplify(T_end_effector[1, 3])
z_position = sp.simplify(T_end_effector[2, 3])

print("X Position:")
print(x_position)
print("Y Position:")
print(y_position)
print("Z Position:")
print(z_position)




# X Position:
# -l2*cos(l1) - l5*(sin(l1)*sin(l3)*sin(zeta_1) - cos(l1)*cos(l3))*cos(l4) + l5*(sin(l1)*sin(zeta_1)*sin(zeta_2)*cos(l3) 
# - sin(l1)*cos(zeta_1)*cos(zeta_2) + sin(l3)*sin(zeta_2)*cos(l1))*sin(l4) + pi*((sin(l1)*sin(l3)*sin(zeta_1) - cos(l1)*cos(l3))*sin(l4) 
# + (sin(l1)*sin(zeta_1)*sin(zeta_2)*cos(l3) 
# - sin(l1)*cos(zeta_1)*cos(zeta_2) + sin(l3)*sin(zeta_2)*cos(l1))*cos(l4))/2 - pi*sin(l1)*sin(zeta_1)/2


# Y Position:
# -l5*(sin(zeta_1)*cos(zeta_2) + sin(zeta_2)*cos(l3)*cos(zeta_1))*sin(l4) + l5*sin(l3)*cos(l4)*cos(zeta_1) 
# - pi*((sin(zeta_1)*cos(zeta_2) + sin(zeta_2)*cos(l3)*cos(zeta_1))*cos(l4) + sin(l3)*sin(l4)*cos(zeta_1))/2 + pi*cos(zeta_1)/2


# Z Position:
# -l2*sin(l1) + l5*(sin(l1)*cos(l3) + sin(l3)*sin(zeta_1)*cos(l1))*cos(l4) + l5*(sin(l1)*sin(l3)*sin(zeta_2) 
# - sin(zeta_1)*sin(zeta_2)*cos(l1)*cos(l3) + cos(l1)*cos(zeta_1)*cos(zeta_2))*sin(l4) 
# - pi*((sin(l1)*cos(l3) + sin(l3)*sin(zeta_1)*cos(l1))*sin(l4) - (sin(l1)*sin(l3)*sin(zeta_2) 
# - sin(zeta_1)*sin(zeta_2)*cos(l1)*cos(l3) + cos(l1)*cos(zeta_1)*cos(zeta_2))*cos(l4))/2 + pi*sin(zeta_1)*cos(l1)/2 + pi/2