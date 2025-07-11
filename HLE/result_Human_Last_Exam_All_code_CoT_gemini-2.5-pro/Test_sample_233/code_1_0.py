def solve_surface_genus():
    """
    This function calculates the smallest positive integer g for the described topology problem.
    """
    
    # Step 1 & 2: Define the relationship between the genera.
    # Let h be the genus of the initial surface Σ.
    # Let h_cap be the genus of the surface S_cap used to cap the boundary of Σ.
    # Let g be the genus of the resulting closed surface Σ' = Σ ∪ S_cap.
    # By gluing the surfaces along their boundary, their genera add up.
    # The formula is: g = h + h_cap.
    
    # The problem states the genus of Σ is 10.
    h = 10
    
    # Step 3 & 4: Determine the maximum required genus for the capping surface.
    # The problem asks for the smallest integer g such that for ANY choice of Σ,
    # a containing surface Σ' of genus g exists.
    # This means we must prepare for the "worst-case" embedding of Σ,
    # which is the one that forces the minimal capping genus, h_cap, to be as large as possible.

    # A theorem in topology states that the minimal genus of a capping surface (a relative Seifert surface)
    # is bounded by the genus of the original surface.
    # So, h_cap <= h. In our case, h_cap <= 10.

    # It is also possible to construct a surface Σ of genus 10 that requires a capping surface of genus 10.
    # This happens when the 10 handles of Σ are embedded in a way that they all pierce through any simple
    # disk that could cap the boundary. To create a valid non-intersecting cap, one must add a handle
    # to the capping surface for each piercing handle of Σ.
    
    # Therefore, the maximum value the minimal h_cap can be forced to take is 10.
    h_cap_max = 10
    
    # Step 5: Calculate the final genus g.
    g = h + h_cap_max
    
    # Step 6: Print the result and the equation.
    print("The problem is to find the smallest genus 'g' of a closed surface Σ' that can contain any given surface Σ.")
    print("Σ has genus h=10 and one unknotted boundary.")
    print("The final genus 'g' is the sum of the genus of Σ and the genus of the capping surface, h_cap.")
    print(f"The genus of the given surface Σ is h = {h}.")
    print(f"The maximum necessary genus for the capping surface is h_cap = {h_cap_max}.")
    print("The final equation is g = h + h_cap.")
    print(f"Thus, the smallest integer g is {h} + {h_cap_max} = {g}.")

solve_surface_genus()