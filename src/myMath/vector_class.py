import numpy as np


class Vector:
    def __init__(self, x, y, z):
        self.x = x
        self.y = y
        self.z = z
        
        
    def __str__(self):
        """Human-readable string representation of the vector."""
        return '({:g})ex + ({:g})ey + ({:g})ez'.format(self.x, self.y, self.z)
    
    def __add__(self, vector):
        x = self.x + vector.x
        y = self.y + vector.y
        z = self.z + vector.z
        return Vector(x, y, z)
       
    def __sub__(self, vector):
        x = self.x - vector.x
        y = self.y - vector.y
        z = self.z - vector.z
        return Vector(x, y, z)
    
    def __mul__(self, scalar):
        x = scalar * self.x
        y = scalar * self.y
        z = scalar * self.z
        return Vector(x, y , z)
    
    def __rmul__(self, scalar):
        x = scalar * self.x
        y = scalar * self.y
        z = scalar * self.z
        return Vector(x, y , z)

    def __truediv__(self, scalar):
        x = self.x / scalar
        y = self.y / scalar
        z = self.z / scalar
        return Vector(x, y, z)
    
    def __pos__(self):
        return self
       
    def __neg__(self):
        return Vector(-self.x, -self.y, -self.z)
    
    # def __rmatmul__(self, A):
    #     x = A[0][0] * self.x + A[0][1] * self.y + A[0][2] * self.z
    #     y = A[1][0] * self.x + A[1][1] * self.y + A[1][2] * self.z
    #     z = A[2][0] * self.x + A[2][1] * self.y + A[2][2] * self.z
    #     return Vector(x, y, z)
    
    def transform(self, A):
        x = A[0][0] * self.x + A[0][1] * self.y + A[0][2] * self.z
        y = A[1][0] * self.x + A[1][1] * self.y + A[1][2] * self.z
        z = A[2][0] * self.x + A[2][1] * self.y + A[2][2] * self.z
        return Vector(x, y, z)
    
    def dot(self, vector):
        return self.x * vector.x + self.y * vector.y + self.z * vector.z
    
    def cross(self, vector):
        x = self.y * vector.z - self.z * vector.y
        y = - (self.x * vector.z - self.z * vector.x)
        z = self.x * vector.y - self.y * vector.x
        return Vector(x, y , z)
    
    def norm(self):
        return np.sqrt(self.dot(self))


if __name__=='__main__':
    
    
    pass