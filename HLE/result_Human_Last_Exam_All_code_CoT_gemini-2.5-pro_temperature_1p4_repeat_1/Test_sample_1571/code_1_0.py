def calculate_damage():
    """
    Calculates the maximum damage Player A can deal to Player B in one turn.
    """
    print("Player A's optimal plan to maximize damage:\n")

    # Step 1: Main Phase Action
    etb_damage = 2
    print(f"1. Player A casts Twinshot Sniper. Its ability triggers upon entering the battlefield.")
    print(f"   - Damage from Twinshot Sniper's ability: {etb_damage}\n")

    # Step 2: Declare Attackers
    # Player A crews Mukotai Soulripper by tapping the two 1/1 Ninja tokens.
    attackers = {
        "Replication Specialist": {"power": 2, "blocked": False},
        "Scrap Welder": {"power": 3, "blocked": False},
        "Ironhoof Boar": {"power": 4, "blocked": False, "must_block": True},
        "Twinshot Sniper": {"power": 2, "blocked": False},
        "Mukotai Soulripper": {"power": 4, "blocked": False, "menace": True}
    }
    print("2. Player A declares attacks with all possible creatures:")
    for name, stats in attackers.items():
        print(f"   - {name} (Power: {stats['power']})")
    print("\n")

    # Step 3: Player B's Optimal Blocking
    print("3. Player B must block to minimize damage and has two 1/1 blockers.")
    # Handle mandatory block
    blockers_available = 2
    for name, stats in attackers.items():
        if stats.get("must_block"):
            stats["blocked"] = True
            blockers_available -= 1
            print(f"   - Due to the 'must be blocked' rule, Player B is forced to use one blocker on Ironhoof Boar, preventing {stats['power']} damage.")

    # Handle remaining blocks
    # Player B blocks the highest power creature they can with their remaining blocker.
    # They cannot block the creature with Menace.
    remaining_attackers_to_block = {name: stats for name, stats in attackers.items() if not stats["blocked"] and not stats.get("menace")}
    if blockers_available > 0 and remaining_attackers_to_block:
      # Find the creature with the highest power to block next
      best_block_target = max(remaining_attackers_to_block, key=lambda x: remaining_attackers_to_block[x]['power'])
      attackers[best_block_target]["blocked"] = True
      blockers_available -= 1
      print(f"   - With their one remaining blocker, Player B blocks the highest-power creature they can: {best_block_target}, preventing {attackers[best_block_target]['power']} damage.\n")

    # Step 4: Calculate Combat Damage
    print("4. The remaining unblocked creatures deal combat damage:")
    combat_damage = 0
    unblocked_creatures = []
    for name, stats in attackers.items():
        if not stats["blocked"]:
            damage = stats['power']
            combat_damage += damage
            unblocked_creatures.append({"name": name, "damage": damage})
            print(f"   - {name} is unblocked and deals {damage} damage.")

    total_damage = etb_damage + combat_damage
    print("\n")

    # Step 5: Final Calculation
    print("5. The total life Player B loses is the sum of the initial ability damage and the combat damage.")
    damage_sources = [etb_damage] + [c['damage'] for c in unblocked_creatures]
    equation = " + ".join(map(str, damage_sources))
    print(f"Final equation: {etb_damage} (Twinshot Sniper ETB) + {unblocked_creatures[0]['damage']} ({unblocked_creatures[0]['name']}) + {unblocked_creatures[1]['damage']} ({unblocked_creatures[1]['name']}) + {unblocked_creatures[2]['damage']} ({unblocked_creatures[2]['name']}) = {total_damage}")

    print(f"\nTotal life lost by Player B: {total_damage}")


calculate_damage()