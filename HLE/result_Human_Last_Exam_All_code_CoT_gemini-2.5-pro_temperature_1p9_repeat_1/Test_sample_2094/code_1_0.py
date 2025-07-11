def solve_hook_flash_question():
    """
    Analyzes which champion(s) can perform a hook-flash, defined as hooking in
    one direction and flashing in the opposite direction to increase range.

    - Blitzcrank, Nautilus, Pyke: Flashing changes the hook's origin. Flashing
      backward while casting forward would *decrease* the hook's effective range.
    - Thresh: Due to his Q's wind-up mechanic, he can cast Q, Flash backward,
      and the hook will travel from his new location to the original target point,
      thereby *increasing* its range.
    """
    champions = ["Thresh"]
    answer = ",".join(champions)
    print(answer)

solve_hook_flash_question()