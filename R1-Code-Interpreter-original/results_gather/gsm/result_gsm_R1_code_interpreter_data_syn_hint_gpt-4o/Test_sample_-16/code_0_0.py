# Number of sisters and their age
num_sisters = 2620008
age_of_sisters = 16

# Age of the brothers
age_of_younger_brother = 12
age_of_older_brother = age_of_younger_brother * 2

# Total age of all sisters
total_age_sisters = num_sisters * age_of_sisters

# Total age of all siblings
total_age_siblings = total_age_sisters + age_of_younger_brother + age_of_older_brother

print(total_age_siblings)