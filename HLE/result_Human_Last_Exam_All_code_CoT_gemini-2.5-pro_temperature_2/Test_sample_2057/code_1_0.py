import collections

# The proposed distribution of hats
# BWBWBWBBW: 5 Black, 4 White
hat_config = ['B', 'W', 'B', 'W', 'B', 'W', 'B', 'B', 'W']
N = len(hat_config)
TOTAL_B = hat_config.count('B')
TOTAL_W = hat_config.count('W')

def get_view_and_unseen(config, person_idx):
    """
    Calculates the number of observed hats and unseen hats for a given person.
    """
    # Person can't see self, left neighbor, right neighbor
    unseen_indices = {(person_idx - 1 + N) % N, person_idx, (person_idx + 1) % N}
    
    seen_b = 0
    seen_w = 0
    for i, hat in enumerate(config):
        if i not in unseen_indices:
            if hat == 'B':
                seen_b += 1
            else:
                seen_w += 1
                
    unseen_b = TOTAL_B - seen_b
    unseen_w = TOTAL_W - seen_w
    return (seen_b, seen_w), (unseen_b, unseen_w)

def check_r1_speaker(config):
    """
    Checks if anyone would speak in Round 1 for a given config.
    Returns the index of a speaker, or -1 if no one speaks.
    """
    for i in range(N):
        seen, unseen = get_view_and_unseen(config, i)
        # If all 3 unseen hats must be of one color, the person knows their color
        if unseen[0] == 3 and unseen[1] == 0:
            return i # This person knows they are Black
        if unseen[0] == 0 and unseen[1] == 3:
            return i # This person knows they are White
    return -1

def is_r1_valid(config):
    """ A configuration is R1-valid if nobody speaks in Round 1. """
    return check_r1_speaker(config) == -1
    
def check_r2_speaker(config):
    """
    Checks if anyone would speak in Round 2 for a given config.
    Returns the index of a speaker, or -1 if no one speaks.
    """
    if not is_r1_valid(config):
        return -1 # Premise is that config is R1-valid
        
    for i in range(N):
        my_hat = config[i]
        other_hat = 'W' if my_hat == 'B' else 'B'
        
        # Consider the hypothesis that the hat is the "other" color
        # This requires reconstructing possible worlds.
        seen, unseen = get_view_and_unseen(config, i)
        
        # If I hypothesize my hat is 'other_hat', my neighbors must have the
        # remaining unseen hats.
        
        # Calculate unseen hats for the hypothetical
        hypo_unseen_neighbors = list(collections.Counter(unseen) - collections.Counter(other_hat))
        
        # There are two ways my neighbors could have these hats (L,R) or (R,L)
        # We need to check if ALL possibilities lead to R1-invalid worlds.
        
        valid_hypo_worlds = 0
        neighbor_perms = [hypo_unseen_neighbors]
        if hypo_unseen_neighbors[0] != hypo_unseen_neighbors[1]:
             neighbor_perms.append(hypo_unseen_neighbors[::-1])

        for perm in neighbor_perms:
            hypo_config = list(config)
            hypo_config[i] = other_hat
            hypo_config[(i - 1 + N) % N] = perm[0]
            hypo_config[(i + 1) % N] = perm[1]
            if is_r1_valid(hypo_config):
                valid_hypo_worlds += 1
                break # Found one valid world for this hypothesis, so can't decide
        
        # If there are ZERO valid R1 worlds for the counter-factual hypothesis, I know my color.
        if valid_hypo_worlds == 0:
            return i # This person would speak in R2
    
    return -1

def is_r2_valid(config):
    return check_r2_speaker(config) == -1
    
def find_r3_speakers(config):
    """
    Finds speakers in Round 3.
    They speak if one of their hypos creates a world that is NOT R2-valid.
    """
    speakers = []
    for i in range(N):
        my_hat = config[i]
        
        # Let's test the counter-factual hypothesis
        other_hat = 'W' if my_hat == 'B' else 'B'
        
        # Build the hypothetical world based on my view and the "other" hat color
        seen, unseen = get_view_and_unseen(config, i)
        hypo_unseen_self = collections.Counter(other_hat)
        
        # Check if this hypo is even possible based on unseen hats
        if (hypo_unseen_self - collections.Counter(unseen)): # Is other_hat in unseen pool?
           continue # Cannot form this hypothesis
           
        hypo_unseen_neighbors_counts = collections.Counter(unseen) - hypo_unseen_self
        hypo_unseen_neighbors = list(hypo_unseen_neighbors_counts.elements())
        
        
        can_rule_out_hypo = True
        
        # Generate possible configurations from the counter-factual hypothesis
        neighbor_perms = [hypo_unseen_neighbors]
        if len(hypo_unseen_neighbors) > 1 and hypo_unseen_neighbors[0] != hypo_unseen_neighbors[1]:
             neighbor_perms.append(hypo_unseen_neighbors[::-1])
        
        for perm in neighbor_perms:
            hypo_config = list(config)
            hypo_config[i] = other_hat
            hypo_config[(i - 1 + N) % N] = perm[0]
            hypo_config[(i + 1) % N] = perm[1]
            
            # This is a potential world under my hypothesis. Is it R1-valid?
            if not is_r1_valid(hypo_config):
                continue # This world would have been disproven in R1
            
            # Is this potential world R2-valid?
            if is_r2_valid(hypo_config):
                # If ANY possible world under my hypo is R2-valid, I cannot rule it out.
                can_rule_out_hypo = False
                break
        
        if can_rule_out_hypo:
            speakers.append(i)
            
    return speakers

print("Proposed hat configuration: B W B W B W B B W\n")

# Verify everyone is silent in R1 and R2 for the initial configuration
is_valid_r1 = is_r1_valid(hat_config)
print(f"Is configuration R1-valid (no one speaks)? {is_valid_r1}")

is_valid_r2 = is_r2_valid(hat_config)
print(f"Is configuration R2-valid (no one speaks)? {is_valid_r2}")
print("-" * 20)

if is_valid_r1 and is_valid_r2:
    r3_speakers = find_r3_speakers(hat_config)
    num_speakers = len(r3_speakers)
    print(f"After three rounds, we need to find who speaks up.")
    print(f"The number of people who replied 'Yes' is: {num_speakers}")
    print(f"\nThe final distribution of hats around the table is:")
    final_equation_str = ' '.join(hat_config)
    print(final_equation_str)
else:
    print("The proposed configuration is not valid for the first two rounds.")
