---
layout: post
title: "Particle based approach for simulating flow over a cylinder"
subtitle: ""
date: 2020-05-29 10:45:13 -0400
background: '/img/b.png'
thumbnail: '/img/icon.png'
---
 In the mainstream CFD, mostly, mesh-based methods are used. Such methods are widely used in the engineering community. But these methods have some drawbacks as generating meshes is a bit cumbersome work and mainly affects the performance of your simulation results. In this regard, meshless approaches are more convenient. 

There are multiple ways to simulate using particle methods, and I will be talking about the Vortex Blob method for 2D problems. In which, we have particles (or blobs) that have vorticity, and they interact in a way so as to simulate the real fluid flow.
## Getting Started
Everything starts with our good old friend Navier Stokes Equation. Just to refresh your memory it is..

$$
\rho \frac{\partial\bar{v}}{\partial t} + \rho (\bar{v}\cdot\nabla) \bar{v} = -\nabla p + \mu \nabla^2 \bar{v}
$$

We have neglected gravity here. And as we are using the vortex blob method, we need to represent this in the form of vorticity and velocity. After taking the curl of the  NS equation we get-

$$
\underbrace{\frac{D\bar{\omega}}{Dt}}_{advection} =\underbrace{\bar{\omega}\nabla{\bar{v}}}_{stretching} - \underbrace{\nu\nabla^2 \bar{\omega}}_{diffusion}
$$

As mentioned below each term, the first term is an advection term second one is for vorticity stretching, and the third one is for diffusion. Also, for 2D problems, vorticity stretching is zero. Now for simplicity, if we consider an inviscid case, then we have to solve  $$ \frac{D\bar{\omega}}{Dt} =0 $$ on the given domain.Now, if we discretize the vorticity on the domain using particles, then we can simulate inviscid flow. Still, there is a problem with how do we find the velocity of each particle given the vorticity field. Without getting into the mathematics, I would directly state here.

$$
(u, v) = k_\delta * \omega
$$

here $$k_\delta$$ is called kernel. It depends on boundary conditions in the flow under consideration. For point vortex with no boundaries kernel is 
$$
\frac{1}{2\pi r^2}(-y, x)
$$

Now for a system of N particles, the velocity of the ith particle is given by representing convolution integral by summation.The integral and summation are mentioned below. Please note that integral over vorticity is replaced with circulation.

$$
(u,v) = \iint k(x-x', y-y') \omega(x', y') dx'dy'
$$

$$

(u, v) = \sum_{j=1}^{\infty}\Gamma _jk(x-x', y-y')

$$

There is also another problem with using point vortices, singularity. Because of this, we get nonphysical numerical solutions. We can use several other non-singular kernels.

Viscosity is similar to the heat equation; it is just the dissipation of energy. It turns out that we can simulate the effect of viscosity by adding gaussian noise to the position of the system with mean 0 and variance $$ \sqrt{2\nu\Delta t}$$.

Now that we know each contributing terms our new governing equation is

$$ \frac{D\bar{\omega}}{Dt} = \nu \nabla^2\omega $$

We can numerically solve this using operator splitting. It means we solve advection and diffusion separately at each time step. Doing this introduces penalty in terms of an error of the order $$ O(\nu\Delta t)$$, which is acceptable.Now let's dig into the implementation part.

We have not yet discussed implementing the no-slip boundary condition. We add blobs at the control point (points over the boundary at which we satisfied no penetration) having total circulation that is equal to the circulation of panels used to satisfy no penetration. These new blobs satisfy the no-slip at control points. The number of blobs can be varied; constraint is just that their total circulation should add up as mentioned before.

## Implementation
For the advection part, we have to calculate velocities due to all other blobs at each blob. For Krasny blob function for velocity field is given below. Velocity function for Krasny blob is-

$$
(u, v)_i=\sum_{j=1, j\ne i}^{\infty}(-y,x)\frac{1}{2\pi (r_j^2+\delta^2)}
$$
```python

def krasny_vel_at_pt(pt, wpos, gamma, delta, exception_index=-1):
    '''
        returns complex velocity at a point 

        pt: location of blob at velocity is calculated
        wpos: positions of all blobs
        gamma, delta: circulation and delta of all blobs

    '''

    v=0+0j
    dz=-wpos+pt
    for k in prange(len(wpos)):
        if k!= exception_index:
            v= v+ (-1j*dz[k].conjugate()*gamma[k] /(2*np.pi*(  (dz[k]* (dz[k].conjugate())) + delta**2 ) ))

    return v.conjugate()

```
Now let's have a look at overall implementation instead of individual functions. As mentioned earlier, to implement using operator splitting, we use the following algorithm. One iteration is-
1. Find velocity on the boundary of the fluid.
2. Introduce blobs at the boundary to satisfy no slip
3. Advect old blobs.
4. Diffuse all blobs.
5. Reflect or delete all the blobs inside cylinder


```python
def one_iter(v_inf, dt, nu, gmax, npan, b_pos, b_gamma, cpt, cpt_b, delta):
    '''
        one iteration for simulation with viscous flow
    '''
    
    gamma_panel, b_gamma_new, b_pos_new = get_p_and_b(b_pos, b_gamma, delta, v_inf, cpt, cpt_b, panels, gmax)
    
    b_pos=advect(b_pos_new ,b_gamma_new ,b_pos, b_gamma, delta, panels, gamma_panel, v_inf, dt, cpt, cpt_b, gmax)

    reflect(b_pos)
    
    b_pos=np.hstack((b_pos, b_pos_new))

    b_gamma=np.hstack((b_gamma, b_gamma_new))
    print(len(b_pos))
    
    b_pos=diffuse(b_pos, nu, dt)
    reflect(b_pos, 1)

    return b_pos, b_gamma
        
```
## Results

The simulation looks like this!!

<!-- [<img src="/img/master.gif" width="800"/>](/img/master.gif) -->

<iframe width="907" height="397" src="https://www.youtube.com/embed/SE9kUbG_cf4" frameborder="0" allow="accelerometer; autoplay; encrypted-media; gyroscope; picture-in-picture" allowfullscreen></iframe>

Before going to get to that final simulation, I performed boundaryless simulations of vortex sheets and patches have a look at some of them.

<iframe width="873" height="491" src="https://www.youtube.com/embed/s2dbctRlhdQ" frameborder="0" allow="accelerometer; autoplay; encrypted-media; gyroscope; picture-in-picture" allowfullscreen></iframe>

For further details, visit my repository [here]( https://github.com/aerorohit/Vortex_method ).

Thank you!
