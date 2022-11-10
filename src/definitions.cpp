#include"definitions.hpp"

namespace bri {

std::size_t range_hash::operator() (range_t const& range) const 
{
    auto hash1 = std::hash<ulint>()(range.first);
    auto hash2 = std::hash<ulint>()(range.second);
    std::size_t seed = 0;
    seed ^= hash1 + 0x9e3779b9 + (seed << 6) + (seed >> 2);
    seed ^= hash2 + 0x9e3779b9 + (seed << 6) + (seed >> 2);
    return seed;
}

std::string get_time(ulint time){

	std::stringstream ss;

	if(time>=3600){

		ulint h = time/3600;
		ulint m = (time%3600)/60;
        ulint s = (time%3600)%60;

		ss  << time << " seconds. ("<< h << "h " << m << "m " << s << "s" << ")";

	}else if (time>=60){

		ulint m = time/60;
		ulint s = time%60;

		ss << time << " seconds. ("<< m << "m " << s << "s" << ")";

	}else{

		ss << time << " seconds.";

	}

	return ss.str();
}

uchar bitsize(ulint x)
{
    if (x == 0) return 1;
    return 64 - __builtin_clzll(x);
}

//parse pizza&chilli patterns header
void header_error() {
	std::cout << "Error: malformed header in patterns file" << std::endl;
	std::cout << "Take a look here for more info on the file format: http://pizzachili.dcc.uchile.cl/experiments.html" << std::endl;
	exit(0);
}

ulint get_number_of_patterns(std::string header){

	ulint start_pos = header.find("number=");
	if (start_pos == std::string::npos or start_pos+7>=header.size())
		header_error();

	start_pos += 7;

	ulint end_pos = header.substr(start_pos).find(" ");
	if (end_pos == std::string::npos)
		header_error();

	ulint n = std::atoi(header.substr(start_pos).substr(0,end_pos).c_str());

	return n;
}

ulint get_patterns_length(std::string header){

	ulint start_pos = header.find("length=");
	if (start_pos == std::string::npos or start_pos+7>=header.size())
		header_error();

	start_pos += 7;

	ulint end_pos = header.substr(start_pos).find(" ");
	if (end_pos == std::string::npos)
		header_error();

	ulint n = std::atoi(header.substr(start_pos).substr(0,end_pos).c_str());

	return n;

}

}; // namespace bri