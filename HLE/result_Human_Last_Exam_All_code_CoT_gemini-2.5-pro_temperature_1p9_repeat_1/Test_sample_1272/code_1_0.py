import sys
import io

# Set stdout to use UTF-8 encoding to handle Vietnamese characters correctly
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8')

def find_mu_sheng_defeat_location():
  """
  This function identifies and prints the modern Vietnamese province where
  Ming General Mu Sheng experienced his first major defeat during the Lam Sơn Uprising.

  Historical Context:
  In late 1427, during the Lam Sơn Uprising, the Ming Dynasty sent two armies to relieve
  their besieged forces in Vietnam. One was led by Liu Sheng, and the other by Mu Sheng.

  The Battle:
  After Liu Sheng's army was famously destroyed at Chi Lăng Pass, the Vietnamese forces
  led by generals Phạm Văn Xảo and Trịnh Khả moved to intercept Mu Sheng's army, which
  was advancing from Yunnan. They set up an ambush at the Cần Trạm pass (ải Cần Trạm),
  also recorded as Khâu Cấp. Mu Sheng's forces were defeated decisively and forced to flee
  back to Yunnan.

  Location Mapping:
  The historical sites of the Cần Trạm pass and the associated battlefields of Lãnh Câu and Đan Xá
  are located in the territory of what is now Yên Bái Province in modern Vietnam.
  """
  province = "Yên Bái"
  print(f"The current Vietnamese province where Ming General Mu Sheng experienced his first major defeat in the 1427 campaign is: {province}")

find_mu_sheng_defeat_location()