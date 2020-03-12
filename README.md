# SwarmGatheringHalfPlanes

Code for the research paper ''Probabilistic Gathering of Agents with Simple Sensors''.

Each folder represents a different scheme of evolution:
- Codes in the folders discrete and continous/memoryless are those described in the paper.
- Code continuous/memory_stopped corresponds to a rule of motion not described in the paper for which we further assume that agents remain static for the rest of the entire interval if at some point in the interval they sense another agent. This is done by adding a memory bit to each agent. Proof of gathering is almost identical to the one presented in the paper with very minor updates.

For display purposes and early stopping, this code uses a function for computing the smallest enclosing circle created by user FSta, based on an example by Yazan Ahed (yash78@gmail.com), who based his code on a Java applet by Shripad Thite (http://heyoka.cs.uiuc.edu/~thite/mincircle/). The link is unfortunately not accessible anymore. Credit for this function is mainly due to the author of the applet.

Contact address: thomas.dages@cs.technion.ac.il
