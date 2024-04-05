# This code is adapted for our purposes from original code online:
# -----------------------------------------------------------------------------
# From Pytnon to Numpy
# Copyright (2017) Nicolas P. Rougier - BSD license
# More information at https://github.com/rougier/numpy-book
# -----------------------------------------------------------------------------
#!/usr/bin/env python
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.path import Path
from matplotlib.animation import FuncAnimation
from matplotlib.collections import PathCollection
import sys
import os
import random
import math

class MarkerCollection:
    """
    Marker collection
    """

    def __init__(self, n=100):
        v = np.array([(-0.25, -0.25), (+0.0, +0.5), (+0.25, -0.25), (0, 0)])
        c = np.array([Path.MOVETO, Path.LINETO, Path.LINETO, Path.CLOSEPOLY])
        self._base_vertices = np.tile(v.reshape(-1), n).reshape(n, len(v), 2)
        self._vertices = np.tile(v.reshape(-1), n).reshape(n, len(v), 2)
        self._codes = np.tile(c.reshape(-1), n)

        self._scale = np.ones(n)
        self._translate = np.zeros((n, 2))
        self._rotate = np.zeros(n)

        self._path = Path(vertices=self._vertices.reshape(n*len(v), 2),
                          codes=self._codes)
        self._collection = PathCollection([self._path], linewidth=0.5,
                                          facecolor="k", edgecolor="w")

    def update(self):
        n = len(self._base_vertices)
        self._vertices[...] = self._base_vertices * self._scale
        cos_rotate, sin_rotate = np.cos(self._rotate), np.sin(self._rotate)
        R = np.empty((n, 2, 2))
        R[:, 0, 0] = cos_rotate
        R[:, 1, 0] = sin_rotate
        R[:, 0, 1] = -sin_rotate
        R[:, 1, 1] = cos_rotate
        self._vertices[...] = np.einsum('ijk,ilk->ijl', self._vertices, R)
        self._vertices += self._translate.reshape(n, 1, 2)


class Flock:
    def __init__(self, count=500, width=500, height=250):
        self.width = width
        self.height = height
        self.min_velocity = 0.5
        self.max_velocity = 2.0
        self.max_acceleration = 0.03
        self.velocity = np.zeros((count, 2), dtype=np.float32)
        self.position = np.zeros((count, 2), dtype=np.float32)

        angle = np.random.uniform(0, 2*np.pi, count)
        self.velocity[:, 0] = np.cos(angle)
        self.velocity[:, 1] = np.sin(angle)
        angle = np.random.uniform(0, 2*np.pi, count)
        radius = min(width, height)/2*np.random.uniform(0, 1, count)
        self.position[:, 0] = width/2 + np.cos(angle)*radius
        self.position[:, 1] = height/2 + np.sin(angle)*radius


    def run(self):

        global distance, counter, sep, ali, coh

        position = self.position
        velocity = self.velocity
        min_velocity = self.min_velocity
        max_velocity = self.max_velocity
        max_acceleration = self.max_acceleration
        n = len(position)

        dx = np.absolute(np.subtract.outer(position[:, 0], position[:, 0]))
        dx = np.minimum(dx, self.width-dx)
        dy = np.absolute(np.subtract.outer(position[:, 1], position[:, 1]))
        dy = np.minimum(dy, self.height-dy)
        distance = np.hypot(dx, dy)

        # Compute common distance masks
        mask_0 = (distance > 0)
        mask_1 = (distance < seprad)
        mask_2 = (distance < alirad)
        mask_3 = (distance < cohrad)
        mask_1 *= mask_0
        mask_2 *= mask_0
        mask_3 *= mask_0
        mask_1_count = np.maximum(mask_1.sum(axis=1), 1)
        mask_2_count = np.maximum(mask_2.sum(axis=1), 1)
        mask_3_count = np.maximum(mask_3.sum(axis=1), 1)

        # Separation
        mask, count = mask_1, mask_1_count
        target = np.dstack((dx, dy))
        target = np.divide(target, distance.reshape(n, n, 1)**2, out=target,
                           where=distance.reshape(n, n, 1) != 0)
        steer = (target*mask.reshape(n, n, 1)).sum(axis=1)/count.reshape(n, 1)
        norm = np.sqrt((steer*steer).sum(axis=1)).reshape(n, 1)
        steer = max_velocity*np.divide(steer, norm, out=steer,
                                       where=norm != 0)
        steer -= velocity

        # Limit acceleration
        norm = np.sqrt((steer*steer).sum(axis=1)).reshape(n, 1)
        steer = np.multiply(steer, max_acceleration/norm, out=steer,
                            where=norm > max_acceleration)

        separation = steer

        # Alignment
        # ---------------------------------------------------------------------
        # Compute target
        mask, count = mask_2, mask_2_count
        target = np.dot(mask, velocity)/count.reshape(n, 1)

        # Compute steering
        norm = np.sqrt((target*target).sum(axis=1)).reshape(n, 1)
        target = max_velocity * np.divide(target, norm, out=target,
                                          where=norm != 0)
        steer = target - velocity

        # Limit acceleration
        norm = np.sqrt((steer*steer).sum(axis=1)).reshape(n, 1)
        steer = np.multiply(steer, max_acceleration/norm, out=steer,
                            where=norm > max_acceleration)
        alignment = steer

        # Cohesion
        # ---------------------------------------------------------------------
        # Compute target
        mask, count = mask_3, mask_3_count
        target = np.dot(mask, position)/count.reshape(n, 1)

        # Compute steering
        desired = target - position
        norm = np.sqrt((desired*desired).sum(axis=1)).reshape(n, 1)
        desired *= max_velocity / norm
        steer = desired - velocity

        # Limit acceleration
        norm = np.sqrt((steer*steer).sum(axis=1)).reshape(n, 1)
        steer = np.multiply(steer, max_acceleration/norm, out=steer,
                            where=norm > max_acceleration)
        cohesion = steer

        # ---------------------------------------------------------------------
        acceleration = sep*separation + ali*alignment + coh*cohesion
        velocity += acceleration

        norm = np.sqrt((velocity*velocity).sum(axis=1)).reshape(n, 1)
        velocity = np.multiply(velocity, max_velocity/norm, out=velocity,
                               where=norm > max_velocity)
        velocity = np.multiply(velocity, min_velocity/norm, out=velocity,
                               where=norm < min_velocity)
        position += velocity

        # Wraparound
        #margin = 10
        #turnFactor = 2
        #sizes = (self.width, self.height)
        #for i in range(n):
        #    for j in range(2):
        #        if position[i, j] < margin:
        #            self.velocity[i, j] += turnFactor
        #        elif position[i, j] > sizes[j]-margin:
        #            self.velocity[i, j] -=  turnFactor
        position += (self.width, self.height)
        position %= (self.width, self.height)
        
        
class Boid(object):
    def __init__(self, index, width, height):
        self.index = index
        self.x = random.random()*width
        self.y = random.random()*height
        self.dx = random.random()*10-5
        self.dy = random.random()*10-5
        self.history = []
        
class FlockSquare(object):
    def __init__(self, numBoids, width, height, avoidFactor = 0.05, centeringFactor = 0.005, matchingFactor = 0.05, visualRange=200):
        self.width = width
        self.height = height
        self.numBoids = numBoids
        self.visualRange = visualRange
        self.avoidFactor = avoidFactor
        self.centeringFactor = centeringFactor
        self.matchingFactor = matchingFactor

        self.boids = [Boid(i, self.width, self.height) for i in range(self.numBoids)]
   
    def distance(self, boid1, boid2):
        return math.sqrt(
        (boid1.x - boid2.x) * (boid1.x - boid2.x) +
          (boid1.y - boid2.y) * (boid1.y - boid2.y)
        )


    def nClosestBoids(self, boid, n):
        sorted_boids = sorted(self.boids, key=lambda x, y: self.distance(boid, x) - self.distance(boid, y))
        return sorted_boids[:n + 1]



    def keepWithinBounds(self, boid):
        # Constrain a boid to within the window. If it gets too close to an edge,
        # nudge it back in and reverse its direction.
        margin = 20
        turnFactor = 1

        if boid.x < margin:
            boid.dx += turnFactor

        if boid.x > width - margin:
            boid.dx -= turnFactor

        if boid.y < margin:
            boid.dy += turnFactor

        if boid.y > height - margin:
            boid.dy -= turnFactor
 


    def flyTowardsCenter(self, boid):
        # Find the center of mass of the other boids and adjust velocity slightly to
        # point towards the center of mass.

        centerX = 0
        centerY = 0
        numNeighbors = 0

        for otherBoid in self.boids:
            if self.distance(boid, otherBoid) < self.visualRange:
                centerX += otherBoid.x
                centerY += otherBoid.y
                numNeighbors += 1

        if numNeighbors>0:
            centerX = centerX / numNeighbors
            centerY = centerY / numNeighbors

        boid.dx += (centerX - boid.x) * self.centeringFactor
        boid.dy += (centerY - boid.y) * self.centeringFactor


    def avoidOthers(self, boid):
        # Move away from other boids that are too close to avoid colliding
        minDistance = 20 # The distance to stay away from other boids
        moveX = 0
        moveY = 0
        for otherBoid in self.boids:
            if otherBoid.index != boid.index:
                if self.distance(boid, otherBoid) < minDistance:
                    moveX += boid.x - otherBoid.x
                    moveY += boid.y - otherBoid.y


        boid.dx += moveX * self.avoidFactor
        boid.dy += moveY * self.avoidFactor


    def matchVelocity(self, boid):
        # Find the average velocity (speed and direction) of the other boids and adjust velocity slightly to match.
        

        avgDX = 0
        avgDY = 0
        numNeighbors = 0

        for otherBoid in self.boids:
            if self.distance(boid, otherBoid) < self.visualRange:
                avgDX += otherBoid.dx
                avgDY += otherBoid.dy
                numNeighbors += 1

        if numNeighbors>0:
            avgDX = avgDX / numNeighbors
            avgDY = avgDY / numNeighbors

        boid.dx += (avgDX - boid.dx) * self.matchingFactor
        boid.dy += (avgDY - boid.dy) * self.matchingFactor



    def limitSpeed(self, boid):
        # Speed will naturally vary in flocking behavior, but real animals can't go arbitrarily fast.
        speedLimit = 15

        speed = math.sqrt(boid.dx * boid.dx + boid.dy * boid.dy)
        if speed > speedLimit:
            boid.dx = (boid.dx / speed) * speedLimit
            boid.dy = (boid.dy / speed) * speedLimit



    def run(self):
        for boid in self.boids:
            # Update the velocities according to each rule
            self.flyTowardsCenter(boid);
            self.avoidOthers(boid);
            self.matchVelocity(boid);
            self.limitSpeed(boid);
            self.keepWithinBounds(boid);

            # Update the position based on the current velocity
            boid.x += boid.dx;
            boid.y += boid.dy;
            boid.history.append([boid.x, boid.y])




def updateSquare(*args):
    global flock, collection, trace, filename, distance, counter

    position = np.asarray([[b.x, b.y] for b in flock.boids])
            
    # Flock updating
    flock.run()
    collection._scale = 10
    collection._translate = position
    collection._rotate = -np.pi/2 + np.arctan2(np.asarray([b.dy for b in flock.boids]),
                                               np.asarray([b.dx for b in flock.boids]))
    collection.update()

    #collect positions
    f = open(filename, "a")
    np.savetxt(f, position, fmt='%.8f', newline='\n')
    f.close()



    # Trace updating
    if trace is not None:
        P = position.astype(int)
        trace[height-1-P[:, 1], P[:, 0]] = .75
        trace *= .99
        im.set_array(trace)
            

def update(*args):
    global flock, collection, trace, filename, distance, counter

    # Flock updating
    flock.run()
    collection._scale = 10
    collection._translate = flock.position
    collection._rotate = -np.pi/2 + np.arctan2(flock.velocity[:, 1],
                                               flock.velocity[:, 0])
    collection.update()

    #collect positions
    f = open(filename, "a")
    np.savetxt(f, flock.position, fmt='%.8f', newline='\n')
    f.close()



    # Trace updating
    if trace is not None:
        P = flock.position.astype(int)
        trace[height-1-P[:, 1], P[:, 0]] = .75
        trace *= .99
        im.set_array(trace)



# -----------------------------------------------------------------------------
if __name__ == '__main__':

    global sep, ali, coh, seprad, alirad, cohrad, filename


    if len(sys.argv) == 6:
        n = int(sys.argv[1])
        sep = float(sys.argv[2])
        ali = float(sys.argv[3])
        coh = float(sys.argv[4])
        seprad = 50
        alirad = 100
        cohrad = 150
        filename = str(sys.argv[5])
    elif len(sys.argv) == 9:
        n = int(sys.argv[1])
        sep = float(sys.argv[2])
        seprad = float(sys.argv[3])
        ali = float(sys.argv[4])
        alirad = float(sys.argv[5])
        coh = float(sys.argv[6])
        cohrad = float(sys.argv[7])
        filename = str(sys.argv[8])
    else:
        print('Please input either 5 or 8 arguments: (to allow different input types, alter code at top of main)')
        print('number of points, separation, alignment, cohesion, outfile')
        print('number of points, separation, sep radius, alignment, ali radius, cohesion, coh radius, outfile')
        exit()

    if os.path.exists(filename):
        os.remove(filename)
    if os.path.exists(filename+"_animation.mp4"):
        os.remove(filename+"_animation.mp4")

    with open(filename, 'w') as outfile:
        outfile.write(str(n) + '\n')
    width, height = 500, 250
    flock = FlockSquare(n, width, height, sep, coh, ali)
    collection = MarkerCollection(n)
    trace = None

    fig = plt.figure(figsize=(10, 10*height/width), facecolor="white")
    ax = fig.add_axes([0.0, 0.0, 1.0, 1.0], aspect=1, frameon=False)
    ax.add_collection(collection._collection)
    ax.set_xlim(0, width)
    ax.set_ylim(0, height)
    ax.set_xticks([])
    ax.set_yticks([])


    animation = FuncAnimation(fig, updateSquare, interval=50, frames=999, repeat=True)
    plt.show()
    animation.save(filename + '_animation.gif', fps=20)
    #writergif = animation.PillowWriter(fps=20)
    #animation.save('filename.gif',writer=writergif)
    print("Animation saved")
