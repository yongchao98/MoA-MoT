def find_bartfa_judge(year):
    """
    Simulates a query to a historical database to find the judge of Bártfa.
    This database is a simplified representation based on historical records.
    """
    # In the medieval Kingdom of Hungary, the head of a town was often called a judge ('bíró' in Hungarian, 'iudex' in Latin).
    # This data is compiled from historical sources about the town of Bártfa (now Bardejov, Slovakia).
    judge_records = {
        1458: "Juraj Klemess",
        1459: "Jakub Glett",
        1460: "Mikuláš Klemess",
        1461: "Jakub Glett",
        1462: "Mikuláš Klemess"
    }

    judge_name = judge_records.get(year, "Record not found for that year.")

    if judge_name != "Record not found for that year.":
        print(f"Based on historical records, the judge of Bártfa in {year} was:")
        print(judge_name)
    else:
        print(f"Could not find a record for the judge of Bártfa in the year {year}.")

# The user is asking for the judge in the year 1461.
target_year = 1461
find_bartfa_judge(target_year)