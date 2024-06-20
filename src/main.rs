use rug::{ops::{Pow, SubFrom}, Complete, Rational};

struct Params {
    bit_count: u32,
    stay_rate: Rational,
    state_match: Rational,
}

fn is_set(num: u64, ind: u32) -> bool {
    num & (1 << ind) != 0
}

fn compute_range(start: u64, dist: u64, params: &Params) -> Vec<Rational> {
    let mut state: u64 = 0;
    let mut res = Vec::new();
    let mut output = start ^ (start >> 1);
    let discrepancies = start.count_ones();
    
    let mut switch_rate: Rational = params.stay_rate.clone();
    switch_rate.sub_from(1);
    let staybyswitch = (&params.stay_rate / switch_rate);
    let switchbystay = staybyswitch.recip_ref().complete();

    let mut mismatch_rate: Rational = params.state_match.clone();
    mismatch_rate.sub_from(1);
    let matchbymismatch = (&params.state_match / &mismatch_rate).complete();
    let mismatchbymatch = matchbymismatch.recip_ref().complete();
    
    let mut current = Rational::from((1, 2));
    current *= (&params.stay_rate).pow(params.bit_count - 1).complete();
    current *= (&params.state_match).pow(params.bit_count - discrepancies).complete();
    current *= (&mismatch_rate).pow(discrepancies).complete();
    
    for step in start+1..start+dist {
        let mut total = current.clone();
        for step_int in 1u64..(1 << params.bit_count) {
            let mut current = current.clone();
            let flip = step_int.trailing_zeros();
            let curbit = is_set(state, flip);
            if flip > 0 {
                current *= if curbit == is_set(state, flip-1) {
                    switchbystay.clone()
                } else {
                    staybyswitch.clone()
                };
            }
            if flip < params.bit_count - 1 {
                current *= if curbit == is_set(state, flip+1) {
                    switchbystay.clone()
                } else {
                    staybyswitch.clone()
                };
            }
            current *= if curbit == is_set(output, flip) {mismatchbymatch.clone()} else {matchbymismatch.clone()};
            total += current;
            state ^= (1 << flip);
        }
        res.push(total);
        let flip = step.trailing_zeros();
        current *= if (is_set(output, flip)) {
            matchbymismatch.clone()
        } else {
            mismatchbymatch.clone()
        };
        output ^= 1 << flip;
        assert!(state == 0);
    }
    res
}

fn main() {
    println!("Hello, world!");
}