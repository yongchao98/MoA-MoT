import math

def calculate_maintenance_fluids():
    """
    Calculates pediatric maintenance fluid needs using the Holliday-Segar method
    and prints the step-by-step calculation.
    """
    weight_kg = 25

    fluid1 = 0
    fluid2 = 0
    fluid3 = 0

    # Fluid for the first 10 kg
    weight1 = min(weight_kg, 10)
    fluid1 = weight1 * 100
    print(f"For the first 10 kg: {weight1} kg * 100 ml/kg = {fluid1} ml")

    # Fluid for the next 10 kg
    if weight_kg > 10:
        weight2 = min(weight_kg - 10, 10)
        fluid2 = weight2 * 50
        print(f"For the next 10 kg: {weight2} kg * 50 ml/kg = {fluid2} ml")

    # Fluid for the remaining weight
    if weight_kg > 20:
        weight3 = weight_kg - 20
        fluid3 = weight3 * 20
        print(f"For the remaining {weight3} kg: {weight3} kg * 20 ml/kg = {fluid3} ml")

    # Total daily fluid
    total_daily_fluid = fluid1 + fluid2 + fluid3
    print(f"\nTotal daily maintenance fluid: {fluid1} ml + {fluid2} ml + {fluid3} ml = {total_daily_fluid} ml")

    # Hourly rate
    hourly_rate = total_daily_fluid / 24.0
    # Use math.ceil to round up to two decimal places to be safe in a clinical context,
    # then format to ensure two decimal places are shown.
    hourly_rate_rounded = math.ceil(hourly_rate * 100) / 100.0
    print(f"Hourly rate: {total_daily_fluid} ml / 24 hours = {hourly_rate_rounded:.2f} cc/hr")

    return hourly_rate_rounded

if __name__ == "__main__":
    final_answer = calculate_maintenance_fluids()
    print(f"\n<<<The calculated maintenance fluid rate is {final_answer:.2f} cc/hr.>>>")
