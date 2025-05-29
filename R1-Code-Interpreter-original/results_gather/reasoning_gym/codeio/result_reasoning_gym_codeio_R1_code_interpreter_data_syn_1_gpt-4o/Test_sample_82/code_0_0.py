import datetime

def is_vehicle_allowed(plate_number, date, time):
    valid_date_format = '%d-%m-%Y'
    valid_time_format = '%I:%M%p'
    
    date_object = datetime.datetime.strptime(date, valid_date_format)
    weekday = date_object.weekday()
    
    # Weekend - no restrictions (5,6 = Saturday,Sunday)
    if weekday == 5 or weekday == 6:
        return True

    time_object = datetime.datetime.strptime(time, valid_time_format)
    time = time_object.time()

    morning_begin_time = datetime.time(6)
    morning_end_time = datetime.time(8, 30)
    evening_begin_time = datetime.time(15)
    evening_end_time = datetime.time(19, 30)

    if (time >= morning_begin_time and time <= morning_end_time) or (time >= evening_begin_time and time <= evening_end_time):
        # Only in peak time check vehicle number
        day = date_object.day
        is_weekday_even = ((day % 2) == 0)
        is_plate_number_even = ((int(plate_number[-1]) % 2) == 0)

        return (is_plate_number_even == is_weekday_even)

    return True

# Test the scenario
plate_number = "abcd 09 8890"
date = "29-04-2022"
time = "9:00PM"

print(is_vehicle_allowed(plate_number, date, time))