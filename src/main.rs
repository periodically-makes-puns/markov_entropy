use std::{cmp::min, collections::BTreeMap};
use histo::Histogram;
use rug::{ops::{CompleteRound, Pow, SubFrom}, Complete, Float, Rational};

struct Params {
    bit_count: u32,
    stay_rate: Rational,
    state_match: Rational,
    buckets: usize
}

fn is_set(num: u64, ind: u32) -> bool {
    num & (1 << ind) != 0
}

fn compute_range(start: u64, dist: u64, params: &Params) -> Vec<(String, Rational)> {
    let mut state: u64 = 0;
    let mut res = Vec::new();
    let mut output = start ^ (start >> 1);
    let discrepancies = start.count_ones();
    
    let mut switch_rate: Rational = params.stay_rate.clone();
    switch_rate.sub_from(1);
    let staybyswitch = &params.stay_rate / switch_rate;
    let switchbystay = staybyswitch.recip_ref().complete();

    let mut mismatch_rate: Rational = params.state_match.clone();
    mismatch_rate.sub_from(1);

    let mut current = Rational::from((1, 2));
    current *= (&params.stay_rate).pow(params.bit_count - 1).complete();
    current *= (&params.state_match).pow(params.bit_count - discrepancies).complete();
    current *= (&mismatch_rate).pow(discrepancies).complete();
    
    let matchbymismatch = &params.state_match / mismatch_rate;
    let mismatchbymatch = matchbymismatch.recip_ref().complete();

    // dbg!(staybyswitch.clone(), switchbystay.clone(), matchbymismatch.clone(), mismatchbymatch.clone());
    // dbg!(current.clone());
    for step in (start+1)..(start+dist+1) {
        let mut total = Rational::ZERO.clone();
        for step_int in 1u64..(1 << params.bit_count)+1 {
            let flip = min(step_int.trailing_zeros(), params.bit_count-1);
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
            state ^= 1 << flip;
            total += current.clone();
        }
        res.push((format!("{:01$b}", output, params.bit_count.try_into().unwrap()), total));
        let flip = min(step.trailing_zeros(), params.bit_count - 1);
        current *= if is_set(output, flip) {
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
    let mut params = Params {
        bit_count: 12,
        stay_rate: Rational::parse("2/3").unwrap().complete(),
        state_match: Rational::parse("1/1000").unwrap().complete(), 
        buckets: 32
    };
    params.state_match.sub_from(1);
    println!("n = {}", params.bit_count);
    println!("P(state unchanged) = {}", params.stay_rate.clone());
    println!("P(bit matches state) = {}", params.state_match.clone());
    let res = compute_range(0, 1 << (params.bit_count), &params);
    let mut entries: Vec<Float> = Vec::new();
    for (_, probability) in res.iter() {
        let res = Float::with_val(1000, probability);
        let lg2 = res.clone().log2();
        entries.push(res * lg2);
    }
    let entropy = Float::sum(entries.iter()).complete(53);
    println!("entropy = {}", -entropy.clone());
    let entropy_rate = (-entropy - 1) / (params.bit_count - 1);
    println!("entropy rate = {}", entropy_rate);
    println!("sum of probabilities = {}", Rational::sum(res.iter().map(|a| &a.1)).complete());
    let mut normalised = Vec::new();
    for i in 0..res.len() {
        if i % 2 == 1 {continue}
        let sp = (&res[i].1 + &res[i+1].1).complete();
        if res[i].1 > res[i+1].1 {
            normalised.push((res[i+1].1.clone() / sp.clone()));
        } else {
            normalised.push((res[i].1.clone() / sp.clone()));
        }
    }
    normalised.sort_unstable();
    // for (prob_lower, total_prob) in normalised {
    //     println!("{} {}", Float::with_val(53, prob_lower), Float::with_val(53, total_prob));
    // }
    println!("bounds = {} to {}", Float::with_val(53, &normalised[0]), Float::with_val(53, &normalised[normalised.len()-1]));
    let lo = normalised[0].clone();
    let hi = normalised[normalised.len()-1].clone();
    let width = (hi - &lo) / params.buckets;
    println!("width = {}", Float::with_val(53, &width));
    let mut cur = lo + &width;
    let mut c = 0;
    // let mut count: Vec<Rational> = Vec::new();
    // count.resize(params.buckets, Rational::ZERO.clone());
    let mut count = Histogram::with_buckets(params.buckets.try_into().unwrap());
    for probability in normalised.iter() {
        while probability >= &cur {
            c += 1;
            cur += &width;
        }
        count.add(c);
    }
    // for i in 0..params.buckets {
    //     println!("{}", count[i]);
    // }
    println!("{}", count);
    let mut appearances: BTreeMap<Rational, usize> = BTreeMap::new();
    for probability in normalised.iter() {
        if let Some(x) = appearances.get_mut(&probability) {
            *x = *x + 1;
        } else {
            appearances.insert(probability.clone(), 1);
        }
    }
    for (k, v) in appearances.iter() {
        println!("{}x of {}", v, k);
    }

}