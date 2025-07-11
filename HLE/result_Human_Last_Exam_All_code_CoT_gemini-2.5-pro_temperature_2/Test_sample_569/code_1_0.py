import sys

def solve_pinochle_score():
    """
    Calculates the total score for a given Pinochle hand scenario.
    """
    # Step 1 & 2: Define the point values for each meld
    my_run_meld = 150
    my_aces_around_meld = 100
    partner_pinochle_meld = 40
    partner_dix_meld = 2 * 10

    # Step 3: Calculate total meld points
    my_total_meld = my_run_meld + my_aces_around_meld
    partner_total_meld = partner_pinochle_meld + partner_dix_meld
    total_meld_points = my_total_meld + partner_total_meld

    # Step 4: Calculate the total points available from tricks
    # Based on "perfect play" with a dominant hand, all trick points are won.
    # The total value of all counter cards in a 48-card Pinochle deck is 240.
    # Aces (8*11=88), Tens (8*10=80), Kings (8*4=32), Queens (8*3=24), Jacks (8*2=16)
    total_trick_counters = 240
    last_trick_bonus = 10
    total_trick_points = total_trick_counters + last_trick_bonus

    # Step 5: Calculate the final total score
    final_score = total_meld_points + total_trick_points

    # Step 6: Print the detailed breakdown and the final equation
    print("Calculating the final score for the hand...\n")
    print("--- Meld Points ---")
    print(f"My Meld: A run ({my_run_meld}) + 8 aces ({my_aces_around_meld}) = {my_total_meld} points")
    print(f"Partner's Meld: A pinochle ({partner_pinochle_meld}) + two 9s of trump ({partner_dix_meld}) = {partner_total_meld} points")
    print(f"Total Meld: {my_total_meld} + {partner_total_meld} = {total_meld_points} points\n")
    
    print("--- Trick Points ---")
    print("Assuming perfect play, your team wins all tricks.")
    print(f"Points from all counter cards in the deck: {total_trick_counters}")
    print(f"Bonus points for winning the last trick: {last_trick_bonus}")
    print(f"Total Trick Points: {total_trick_counters} + {last_trick_bonus} = {total_trick_points} points\n")

    print("--- Final Score Calculation ---")
    print("The final score is the sum of all meld and trick points.")
    # The final equation with each component number
    equation_str = (
        f"Final Score = {my_run_meld} (my run) + {my_aces_around_meld} (my aces) + "
        f"{partner_pinochle_meld} (partner's pinochle) + {partner_dix_meld} (partner's 9s) + "
        f"{total_trick_counters} (all counters) + {last_trick_bonus} (last trick) = {final_score}"
    )
    print(equation_str)

    # Required for the final answer block
    sys.stdout.write(f"\n<<<{final_score}>>>")

solve_pinochle_score()