import math

print("Yes, your son can calculate the time. Here is a step-by-step calculation using small integers and fractions.")

print("\nStep 1: Find the height of the skyscraper (h).")
print("The angle α is a quarter of a right angle (90 degrees), so α = 90 / 4 = 22.5 degrees.")
print("To find the height, we need the tangent of that angle, tan(22.5°). The exact value is a bit tricky, it's sqrt(2) - 1.")
print("A good and simple fraction for sqrt(2) is 17/12. So we can approximate tan(22.5°) ≈ 17/12 - 1 = 5/12.")
print("The distance (d) is 240 m. To make calculations easier, we can write this as 24 * 10 m.")
print("The height is h = d * tan(α). Let's calculate it:")
print("h ≈ (24 * 10) * (5 / 12)")
print("We can rearrange this to make it simpler: h ≈ (24 / 12) * 10 * 5 = 2 * 5 * 10 m.")

print("\nStep 2: Set up the equation for the fall time (t).")
print("The formula for an object falling is h = (1/2) * g * t^2. We can rewrite this to find the time: t = sqrt(2 * h / g).")
print("For gravity (g), we can use the simple approximation g ≈ 10 m/s^2.")
print("Now let's find the value of t^2 first:")
print("t^2 ≈ (2 * h) / g ≈ (2 * (2 * 5 * 10)) / 10")
print("This simplifies to: t^2 ≈ 2 * 2 * 5 = 20.")

print("\nStep 3: Calculate the final time t.")
print("If t^2 = 20, then t = sqrt(20). This is the same as t = sqrt(4 * 5), which is 2 * sqrt(5).")
print("We need a fraction for sqrt(5). A very close one is 9/4 (which is 2.25).")
print("So, we can find our final answer for t:")
print("t ≈ 2 * (9 / 4)")
print("t ≈ 18 / 4 = 9 / 2 seconds. That's 4.5 seconds.")