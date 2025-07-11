import math

def solve_triangle_problem():
    """
    This function solves for the range of the median m given the height AD and angle bisector AE.
    """
    # Given values
    AD = 12  # Length of the altitude from A to BC
    AE = 13  # Length of the angle bisector from A

    # --- Step 1: Find the length of DE ---
    # Triangle ADE is a right-angled triangle because AD is an altitude (AD ⊥ BC).
    # The right angle is at D. AE is the hypotenuse.
    # By the Pythagorean theorem: AE^2 = AD^2 + DE^2
    de_squared = AE**2 - AD**2
    de = math.sqrt(de_squared)
    
    print("Step 1: Calculate the distance DE.")
    print(f"In the right-angled triangle ADE, we have AD = {AD} and AE = {AE}.")
    print(f"Using the Pythagorean theorem, DE^2 = AE^2 - AD^2.")
    print(f"DE^2 = {AE}^2 - {AD}^2 = {AE**2} - {AD**2} = {de_squared}")
    print(f"So, DE = sqrt({de_squared}) = {de}")
    print("-" * 40)

    # --- Step 2: Express DF in terms of m ---
    # AF is the median of length m. Triangle ADF is also a right-angled triangle (right angle at D).
    # By the Pythagorean theorem: AF^2 = AD^2 + DF^2
    # m^2 = AD^2 + DF^2
    print("Step 2: Establish a relationship for the median m.")
    print("In the right-angled triangle ADF, we have AD = 12 and AF = m.")
    print("Using the Pythagorean theorem, m^2 = AD^2 + DF^2.")
    print(f"m^2 = {AD}^2 + DF^2 = {AD**2} + DF^2")
    print(f"This gives us DF^2 = m^2 - {AD**2}.")
    print("-" * 40)

    # --- Step 3: Apply the condition for angle A to be acute ---
    # A geometric theorem states that angle A is acute if and only if DF > DE.
    # If A were right, DF = DE. If A were obtuse, DF < DE.
    print("Step 3: Use the geometric condition for ∠A to be acute.")
    print("For ∠A to be acute, the distance from the foot of the median to the foot of the altitude (DF)")
    print("must be greater than the distance from the foot of the angle bisector to the foot of the altitude (DE).")
    print("The condition is: DF > DE")
    print("-" * 40)
    
    # --- Step 4: Solve for the range of m ---
    # We have DF > DE, and since both are lengths, we can square both sides: DF^2 > DE^2.
    # Substitute the expressions from Step 1 and 2.
    # m^2 - AD^2 > de_squared
    # m^2 > de_squared + AD^2
    # m^2 > 25 + 144
    # m^2 > 169
    # m > 13
    print("Step 4: Combine the relationships to find the range of m.")
    print(f"From the condition DF > DE, we have DF > {de}.")
    print(f"Squaring both sides gives DF^2 > {de}^2, so DF^2 > {de_squared}.")
    print(f"Substituting DF^2 = m^2 - {AD**2}, we get:")
    print(f"m^2 - {AD**2} > {de_squared}")
    m_squared_rhs = de_squared + AD**2
    print(f"m^2 > {de_squared} + {AD**2}")
    print(f"m^2 > {m_squared_rhs}")
    m_lower_bound = math.sqrt(m_squared_rhs)
    print(f"Taking the square root, we find m > {m_lower_bound}.")
    print("-" * 40)
    print("Final Result:")
    print(f"The range of values for m for which angle A will be acute is m > {m_lower_bound}.")

solve_triangle_problem()