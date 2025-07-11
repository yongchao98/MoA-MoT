def solve_random_walk_speed():
    """
    Solves for the asymptotic speed of a biased random walk on a random ladder graph.

    Problem breakdown:
    1.  The graph is an infinite ladder Z x {0,1} with some edges randomly removed.
        - Lower horizontal edges ((n-1,0),(n,0)) are always present.
        - Upper horizontal edges ((n-1,1),(n,1)) exist with probability 2/3.
        - Vertical edges ((n,0),(n,1)) exist with probability 1/2.
    2.  The random walk transition probability from x to y is proportional to exp(c*(y1-x1)).
    3.  We need to find the limit of the asymptotic speed v(c) as c -> infinity.

    Reasoning for the c -> infinity limit:
    - The bias term exp(c) for forward moves dominates all other transition weights.
    - Let's analyze the vertical movement of the walker.

    - Transition from Level 0 to Level 1:
      A walker at (n,0) can move to (n+1,0) with weight exp(c) or to (n,1) with weight 1 (if the edge exists).
      The probability of moving up is P(up) = 1 / (exp(c) + ...), which tends to 0 as c -> infinity.
      Therefore, the lower level is an "absorbing" state. Once the walker is on level 0, it stays there.

    - Transition from Level 1 to Level 0:
      A walker at (n,1) will move to (n+1,1) if that edge exists. If not (with probability 1/3), it will
      look for other options. The next best option is moving to (n,0) (weight 1), which is infinitely
      preferred over moving backward (weight exp(-c)).
      So, the walker can and will move from level 1 to level 0.

    - Long-term behavior:
      Since the walker can only move from level 1 down to 0 (and not up), any walker will eventually
      end up on level 0 as time t -> infinity. The stationary distribution of the walker is entirely
      concentrated on level 0.

    - Speed Calculation:
      The asymptotic speed is determined by the walker's behavior on level 0.
      On level 0, all horizontal edges exist. The walker at (n,0) moves to (n+1,0) with a
      probability that approaches 1 as c -> infinity.
      Each step takes 1 unit of time and provides 1 unit of horizontal displacement.
      Speed = Displacement / Time = 1 / 1 = 1.
    """
    speed = 1
    # The final equation is simply Speed = 1 / 1
    numerator = 1  # displacement per step on the lower level
    denominator = 1 # time per step on the lower level
    print(f"The asymptotic speed v(c) is determined by the behavior on the absorbing lower level.")
    print(f"On the lower level, the walker moves 1 unit to the right in 1 time step.")
    print(f"Speed = Displacement / Time")
    print(f"v = {numerator} / {denominator}")
    print(f"The limit of the asymptotic speed is: {speed}")

solve_random_walk_speed()