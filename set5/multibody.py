from __future__ import print_function

# multibody.py
# Aritra Biswas
#
# Framework for evolving motion of multiple particles under mutual
# interactions, using the Barnes-Hut tree algorithm.

from Vector import *
import math
import visual as vs
import random
import time
import cPickle as pickle
import copy

class Particle:
    '''A particle in 3D space, holding a position and velocity.'''

    def __init__(self,
            pos = Vector([0.0, 0.0, 0.0]),
            vel = Vector([0.0, 0.0, 0.0]),
            mass = 1.0,
            radius = 0.1,
        ):

        assert isinstance(pos, Vector) and len(pos) == 3
        assert isinstance(vel, Vector) and len(vel) == 3

        self.pos = pos
        self.vel = vel
        self.mass = mass
        self.radius = radius
    
    def set_pos(self, pos):
        '''Set the position in Cartesian coordinates.'''
        
        assert isinstance(pos, Vector) and len(pos) == 3
        self.pos = pos
    
    def set_vel(self, vel):
        '''Set the velocity in Cartesian coordinates.'''

        assert isinstance(vel, Vector) and len(vel) == 3
        self.vel = vel

    def set_pos_sph(self, pos_sph):
        '''Set the position in spherical coordinates, where the
        input pos_ph = [r, theta, phi]. The physics convention is
        used, where phi is the azimuthal angle.'''

        assert isinstance(pos_sph, Vector) and len(pos_sph) == 3
        (r, theta, phi) = list(pos_sph)

        x = r * math.sin(theta) * math.cos(phi)
        y = r * math.sin(theta) * math.sin(phi)
        z = r * math.cos(theta)

        pos = Vector([x, y, z])
        self.set_pos(pos)

    def set_vel_sph(self, vel_sph):
        '''Set the velocity in spherical coordinates, where the
        input vel_sph = [v, theta, phi]. The physics convention is
        used, where phi is the azimuthal angle.'''

        assert isinstance(vel_sph, Vector) and len(vel_sph) == 3
        (v, theta, phi) = list(vel_sph)

        vx = v * math.sin(theta) * math.cos(phi)
        vy = v * math.sin(theta) * math.sin(phi)
        vz = v * math.cos(theta)

        vel = Vector([vx, vy, vz])
        self.set_vel(vel)

def grav_force(test, source):
    '''Returns the force vector of the particle test due to
    the gravitational force from particle source.'''

    assert isinstance(test, Particle)
    assert isinstance(source, Particle)

    if test == source:
        return Vector([0.0, 0.0, 0.0])

    G = 1

    sep = source.pos - test.pos
    force = G * test.mass * source.mass * sep / math.pow(abs(sep), 3)

    return force

def soft_grav_force(test, source):
    '''Similar to grav_force, except the effective distance between two
    particles is always incremented by a small offset, such that the force
    does not blow up when particles are near each other.'''
    
    assert isinstance(test, Particle)
    assert isinstance(source, Particle)
    
    if test == source:
        return Vector([0.0, 0.0, 0.0])

    G = 1
    a = 0.1

    sep = source.pos - test.pos
    force = G * test.mass * source.mass * sep / (math.pow(abs(sep) + a, 3))

    return force

class QuadNode(object):
    '''A node in a Barnes-Hut quadtree, holding its particles, its center of
    mass, and its total mass.'''

    def __init__(self, particles):

        # list of particles in this node
        self.particles = particles

        # dimensions: bottom-left point (begin), and side length
        self.begin = None
        self.side = None
        
        # center of mass and total mass
        self.cm = None
        self.mass = None

        # children
        self.child_bl = None
        self.child_br = None
        self.child_tl = None
        self.child_tr = None

    def valid(self):
        '''Check if the node has valid geometry values.'''

        vals = [self.begin, self.side]
        return all(val is not None for val in vals)

    def calculated(self):
        '''Check if the CM and mass of the node have been calculated.'''

        vals = [self.cm, self.mass]
        return all(val is not None for val in vals)

    def parent(self):
        '''Check if the node has children.'''

        children = [self.child_bl, self.child_br, self.child_tl,
                    self.child_tr]
        return all(child is not None for child in children)

    def split(self):
        '''If the node holds more than one particle, create children
        recursively.'''

        assert self.valid()
        assert not self.parent()
        
        if len(self.particles) > 1:

            # get midpoint coordinates for splitting
            (begin_x, begin_y, _) = self.begin
            mid_x = begin_x + (self.side / 2)
            mid_y = begin_y + (self.side / 2)

            # group particles
            bl_particles = []
            br_particles = []
            tl_particles = []
            tr_particles = []
            for p in self.particles:
                if p.pos[0] < mid_x:
                    if p.pos[1] < mid_y:
                        bl_particles.append(p)
                    else:
                        tl_particles.append(p)
                else:
                    if p.pos[1] < mid_y:
                        br_particles.append(p)
                    else:
                        tr_particles.append(p)

            # bottom-left child
            self.child_bl = QuadNode(bl_particles)
            self.child_bl.begin = Vector([begin_x, begin_y, 0])
            self.child_bl.side = self.side / 2
            self.child_bl.split()
            
            # bottom-right child
            self.child_br = QuadNode(br_particles)
            self.child_br.begin = Vector([mid_x, begin_y, 0])
            self.child_br.side = self.side / 2
            self.child_br.split()
            
            # top-left child
            self.child_tl = QuadNode(tl_particles)
            self.child_tl.begin = Vector([begin_x, mid_y, 0])
            self.child_tl.side = self.side / 2
            self.child_tl.split()
            
            # top-right child
            self.child_tr = QuadNode(tr_particles)
            self.child_tr.begin = Vector([mid_x, mid_y, 0])
            self.child_tr.side = self.side / 2
            self.child_tr.split()

    def calculate(self):
        '''Recursively calculate CM and mass information.'''

        if not self.parent():
            self.mass = sum(p.mass for p in self.particles)
            if self.mass != 0:
                self.cm = sum(p.mass * p.pos for p in self.particles)
                self.cm = self.cm / self.mass
            else:
                self.cm = Vector([0.0, 0.0, 0.0])
        else:
            children = [self.child_bl, self.child_br, self.child_tl,
                        self.child_tr]
            for child in children:
                child.calculate()
            self.mass = sum(c.mass for c in children)
            self.cm = sum(c.mass * c.cm for c in children) / self.mass

class RootNode(QuadNode):
    '''The root node in the quadtree, with additional initialization
    steps.'''

    def __init__(self, particles):

        # initialize regular node
        super(RootNode, self).__init__(particles)

        # compute beginning point
        begin_x = min(p.pos[0] for p in self.particles)
        begin_y = min(p.pos[1] for p in self.particles)
        self.begin = Vector([begin_x, begin_y, 0])

        # compute side length
        side_x = max(p.pos[0] for p in self.particles) - begin_x
        side_y = max(p.pos[1] for p in self.particles) - begin_y
        self.side = max(side_x, side_y)

class QuadTree:
    '''Creates a Barnes-Hut quadtree from a cluster object.'''

    def __init__(self, cluster):

        self.root = RootNode(cluster.particles)
        self.root.split()
        self.root.calculate()

def calc_force(test, quadnode, cluster, theta):
    '''Calculates the force exerted by the particles in the quadnode
    on test, using cluster parameters with opening angle theta.'''

    assert isinstance(test, Particle)
    assert isinstance(quadnode, QuadNode)
    assert isinstance(cluster, Cluster)
    assert theta >= 0

    force = Vector([0.0, 0.0, 0.0])

    if len(quadnode.particles) == 1:
        force = force + cluster.force(test, quadnode.particles[0])

    elif len(quadnode.particles) > 1:
        d = abs(test.pos - quadnode.cm)
        if d != 0:
            ratio = quadnode.side / d
        else:
            ratio = theta + 1
        if ratio < theta:
            effective = Particle(quadnode.cm,
                                 Vector([0.0, 0.0, 0.0]),
                                 quadnode.mass)
            force = force + cluster.force(test, effective)
        else:
            children = [quadnode.child_bl, quadnode.child_br,
                        quadnode.child_tl, quadnode.child_tr]
            for child in children:
                force = force + calc_force(test, child, cluster, theta)

    return force

class Cluster:
    '''A cluster of particles interacting with each other under some
    force.'''

    # INITIALIZATION AND PARTICLE MEMBERSHIP 

    def __init__(self):

        # start cluster without any member particles
        self.particles = []

        # interaction force between pairs of particles
        # force(test, source) should return force vector on test
        self.force = grav_force

        # keep track of evolution timesteps
        self.timestep = 0

        # Opening angle for Barnes-Hut approximation algorithm. 0 means
        # no clustering; higher values cause more clustering and approximation.
        self.theta = 0

    def add_particle(self, particle):
        '''Add a particle to the cluster.'''

        assert isinstance(particle, Particle)
        self.particles.append(particle)

    # EVOLUTION

    def evolve(self, dt):
        '''Evolve by flipping between two types of symplectic Euler
        steps. Output the current step and time.'''

        if self.timestep % 2 == 0:
            self.evolve_symplectic1(dt)
        else:
            self.evolve_symplectic2(dt)
        
        out = "step = {} \t time = {}".format(
                self.timestep,
                self.timestep * dt)
        print(out, end = "\r")
        
        self.timestep += 1

    def evolve_symplectic1(self, dt):
        '''Using the symplectic Euler method, update the positions of
        all particles in the cluster for a timestep dt. Update
        velocities before positions.'''
        
        quadtree = QuadTree(self)

        for test in self.particles:

            # calculate the net acceleration on test
            force = calc_force(test, quadtree.root, self, self.theta)
            accel = force / test.mass

            # use current acceleration to update velocity
            test.set_vel(test.vel + accel * dt)

            # use updated velocity to update position
            test.set_pos(test.pos + test.vel * dt)

        del quadtree
    
    def evolve_symplectic2(self, dt):
        '''Using the symplectic Euler method, update the positions of
        all particles in the cluster for a timestep dt. Update positions
        before velocities.'''
        
        quadtree = QuadTree(self)

        for test in self.particles:
            
            # use current velocity to update position
            test.set_pos(test.pos + test.vel * dt)

            # calculate the net acceleration on test
            force = calc_force(test, quadtree.root, self, self.theta)
            accel = force / test.mass

            # use updated acceleration to update velocity
            test.set_vel(test.vel + accel * dt)
        
        del quadtree

    # VPYTHON VISUALIZATION

    def vs_init(self):
        '''Set up the visualization with the current particles.'''

        # define a list of acceptable colors to use for particles
        colors = [vs.color.red, vs.color.yellow, vs.color.green,
                  vs.color.orange, vs.color.white, vs.color.blue,
                  vs.color.cyan, vs.color.magenta]

        # initialize a dictionary to link particles with visual objects,
        # such that vs_map[particle] = vs_object
        self.vs_map = {}

        for particle in self.particles:

            # create visual object using particle parameters
            vs_obj = vs.sphere()
            vs_obj.radius = particle.radius
            vs_obj.pos = list(particle.pos)
            vs_obj.color = random.choice(colors)

            # link visual object to particle
            self.vs_map[particle] = vs_obj

    def vs_update(self):
        '''Update the visualization with the current positions of
        the particles.'''

        for particle in self.particles:
            vs_obj = self.vs_map[particle]
            vs_obj.pos = list(particle.pos)

    def vs_run(self, dt, tx):
        '''Continuously evolve and update the visual simulation
        with timestep dt and visual speed tx relative to real-time.'''

        while True:
            self.evolve(dt)
            vs.rate(tx / dt)
            self.vs_update()

    # OUTPUT SAVING

    def io_run(self, dt, steps, wrstep, name = "out.pkl"):
        '''Evolve the cluster for given steps with timestep dt and write
        out the cluster every wrstep timesteps.'''

        f_output = open(name, "wb")

        while self.timestep < steps:

            self.evolve(dt)

            if self.timestep % wrstep == 1:

                obj = copy.deepcopy(self)
                pickle.dump(obj, f_output)
                del obj
                
                print("wrte @")

        f_output.close()
