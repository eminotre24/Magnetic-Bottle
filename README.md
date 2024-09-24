# Magnetic Bottle
This is a code which simulates a system of particles confined in a "magnetic bottle", which is a tool used to confine some types of materials, such as plasma, using magnetic fields.

This repository contains a file named "magbottle.py", which is the library that contains all the functions used for the simulation. And the file "main.py" contains an example of a system of 3 particles.

This simulator is supposed to work for a system of n particles (if the system is stable it shoudnt have problems running it), and was made using the algorithm Velocity Verlet, with an extra step, which is an extra calculation of the velocity (an extra half step), as we need the velocity for calculating the different magnetic interaction between the charged particles, and with the magnetic field. 

All the information regarding the simulation is stored in the dataclasses "particle" and "time" (in their respectively variable), after running the function "velocity_verlet". This information is then used in the functions "tray" and "anim", which obviously would requiere the particles to have a list of positions to animate.
