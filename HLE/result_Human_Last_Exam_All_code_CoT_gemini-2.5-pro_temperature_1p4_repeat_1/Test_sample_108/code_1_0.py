import math

def calculate_min_distance():
    """
    Calculates the minimum distance from a robot arm's shoulder to its finger
    based on segment lengths and joint clearance constraints.
    """
    # Segment lengths in cm
    L1 = 40  # Shoulder to Elbow
    L2 = 28  # Elbow to Wrist
    L3 = 15  # Wrist to Hand
    L4 = 10  # Hand to Finger
    
    # Clearance constraint in cm
    C = 3.5

    print("Step 1: Calculate the minimum internal angle for each joint based on the 3.5 cm clearance rule.")
    
    # Elbow joint (between L1 and L2)
    # The fold is limited by the shorter segment adjacent to the joint.
    # sin(alpha_E) = C / min(L1, L2)
    sin_alpha_E = C / L2
    alpha_E = math.asin(sin_alpha_E)
    print(f"Shoulder-Elbow (L1) = {L1} cm, Elbow-Wrist (L2) = {L2} cm")
    print(f"Calculated sin(alpha_Elbow) = {sin_alpha_E:.4f}, so alpha_Elbow = {math.degrees(alpha_E):.2f} degrees.")
    
    # Wrist joint (between L2 and L3)
    # sin(alpha_W) = C / min(L2, L3)
    sin_alpha_W = C / L3
    alpha_W = math.asin(sin_alpha_W)
    print(f"Elbow-Wrist (L2) = {L2} cm, Wrist-Hand (L3) = {L3} cm")
    print(f"Calculated sin(alpha_Wrist) = {sin_alpha_W:.4f}, so alpha_Wrist = {math.degrees(alpha_W):.2f} degrees.")

    # Hand joint (between L3 and L4)
    # sin(alpha_H) = C / min(L3, L4)
    sin_alpha_H = C / L4
    alpha_H = math.asin(sin_alpha_H)
    print(f"Wrist-Hand (L3) = {L3} cm, Hand-Finger (L4) = {L4} cm")
    print(f"Calculated sin(alpha_Hand) = {sin_alpha_H:.4f}, so alpha_Hand = {math.degrees(alpha_H):.2f} degrees.")
    print("-" * 30)

    print("Step 2: Use the Law of Cosines to find the distance from the shoulder to the wrist (SW).")
    # Triangle SEW (Shoulder-Elbow-Wrist)
    # SW^2 = SE^2 + EW^2 - 2 * SE * EW * cos(alpha_E)
    SW_sq = L1**2 + L2**2 - 2 * L1 * L2 * math.cos(alpha_E)
    SW = math.sqrt(SW_sq)
    print(f"Distance Shoulder-Wrist (SW)^2 = {L1}^2 + {L2}^2 - 2*{L1}*{L2}*cos({math.degrees(alpha_E):.2f}) = {SW_sq:.2f}")
    print(f"Distance Shoulder-Wrist (SW) = {SW:.2f} cm.")
    print("-" * 30)

    print("Step 3: Calculate intermediate angles and then the distance from the shoulder to the hand (SH).")
    # In triangle SEW, find angle SWE using the Law of Sines
    # sin(SWE) / L1 = sin(alpha_E) / SW
    sin_angle_SWE = L1 * math.sin(alpha_E) / SW
    angle_SWE = math.asin(sin_angle_SWE)
    
    # To fold tightly, the angle SWH is the difference between SWE and alpha_W
    angle_SWH = angle_SWE - alpha_W
    
    # Triangle SWH (Shoulder-Wrist-Hand)
    # SH^2 = SW^2 + WH^2 - 2 * SW * WH * cos(SWH)
    SH_sq = SW**2 + L3**2 - 2 * SW * L3 * math.cos(angle_SWH)
    SH = math.sqrt(SH_sq)
    print(f"Angle SWE = {math.degrees(angle_SWE):.2f} degrees.")
    print(f"Angle SWH = |SWE - alpha_Wrist| = {math.degrees(abs(angle_SWH)):.2f} degrees.")
    print(f"Distance Shoulder-Hand (SH)^2 = {SW:.2f}^2 + {L3}^2 - 2*{SW:.2f}*{L3}*cos({math.degrees(abs(angle_SWH)):.2f}) = {SH_sq:.2f}")
    print(f"Distance Shoulder-Hand (SH) = {SH:.2f} cm.")
    print("-" * 30)

    print("Step 4: Calculate more intermediate angles and the final distance from the shoulder to the finger (SF).")
    # In triangle SWH, find angle SHW using the Law of Cosines
    # SW^2 = SH^2 + WH^2 - 2 * SH * WH * cos(SHW)
    cos_angle_SHW = (SH**2 + L3**2 - SW**2) / (2 * SH * L3)
    # Handle potential floating point inaccuracies leading to domain errors
    if cos_angle_SHW > 1.0: cos_angle_SHW = 1.0
    if cos_angle_SHW < -1.0: cos_angle_SHW = -1.0
    angle_SHW = math.acos(cos_angle_SHW)

    # To continue the tight fold, the angle SHF is the difference between SHW and alpha_H
    angle_SHF = angle_SHW - alpha_H
    
    # Triangle SHF (Shoulder-Hand-Finger)
    # SF^2 = SH^2 + HF^2 - 2 * SH * HF * cos(SHF)
    SF_sq = SH**2 + L4**2 - 2 * SH * L4 * math.cos(angle_SHF)
    SF = math.sqrt(SF_sq)
    print(f"Angle SHW = {math.degrees(angle_SHW):.2f} degrees.")
    print(f"Angle SHF = |SHW - alpha_Hand| = {math.degrees(abs(angle_SHF)):.2f} degrees.")
    print(f"Final Distance (SF)^2 = {SH:.2f}^2 + {L4}^2 - 2*{SH:.2f}*{L4}*cos({math.degrees(abs(angle_SHF)):.2f}) = {SF_sq:.2f}")
    print(f"Final Minimum Distance (SF) = {SF:.2f} cm.")

calculate_min_distance()