import math

def solve():
    # The problem asks for the limit of the asymptotic speed v(c) as c -> infinity.
    #
    # Let's analyze the behavior of the random walk for very large c.
    # The transition probabilities are proportional to exp(c * delta_x), where delta_x is the change in the horizontal coordinate.
    # For c -> infinity, exp(c) >> 1 >> exp(-c). This makes the walk deterministic and greedy:
    # 1. If a move to the right is possible, the walker will take it with probability approaching 1.
    # 2. If no move to the right is possible, but a vertical move is, the walker will take the vertical one.
    # 3. If neither a right nor a vertical move is possible, the walker must move left.
    #
    # The graph has two "lanes":
    # - The lower lane (y=0) is a complete, infinite line. A walker at (n,0) can always move to (n+1,0). The weight for this move is exp(c). The weight for moving up to (n,1) is 1. For large c, the walker will always choose to move right. Thus, if a walker ever reaches the lower lane, it will stay there, moving right at each step. The speed in this lane is 1.
    #
    # - The upper lane (y=1) is not complete. Its horizontal edges are missing with probability 1/3. This means the upper lane consists of finite segments (with probability 1).
    #
    # Let's consider a walk.
    # If the walk starts on the lower lane, its speed will be 1.
    # If the walk starts on the upper lane, it will move right as long as the path exists. When it hits the end of a segment (no rightward edge), it must find another way. The greedy choice is to move down. If it can't move down at its current position, it will move left until it finds a position from where it can move down. Since every segment in the infinite component must be connected to the lower lane, such a downward path always exists.
    #
    # Crucially, this means any walk will reach the lower lane in a finite number of steps. Once on the lower lane, it will never leave, as moving up (weight 1) is infinitely less likely than moving right (weight exp(c)).
    #
    # The asymptotic speed is defined as lim_{t->inf} E[x_t - x_0] / t.
    # Since the "transient phase" of reaching the lower lane takes a finite amount of time, its contribution to the total displacement gets washed out as t -> infinity. The long-term speed is determined by the behavior in the absorbing lower lane, which is 1.
    #
    # Therefore, the limit of the asymptotic speed v(c) as c -> infinity is 1.
    
    lim_v_c = 1.0
    
    # We are asked to output each number in the final equation. 
    # The final equation is simply the result.
    
    print(int(lim_v_c))

solve()