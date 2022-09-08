// Knuth's model of a fair die using only fair coins
dtmc

const double p;
const double q;

module die

	// local state
	s : [0..2] init 0;
	
	[] s=0 -> 0 : (s'=0) + p : (s'=1) + (1-p) : (s'=2);
	[] s=1 -> 0.8 : (s'=0) + 0.2 : (s'=1) + 0 : (s'=2);
	[] s=2 -> 0 : (s'=0) + q : (s'=1) + (1-q) : (s'=2);

	
endmodule

rewards "coin_flips"
	[] s<2 : 1;
endrewards

