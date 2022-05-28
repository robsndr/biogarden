
pub fn linear_gap_penalty(a: &usize) -> i32 {
    // print!("{}   {}\n", a - 65 , b-65);
    // BLOSUM62[((*a as usize) - 65, (*b as usize) - 65)]
    if *a == 0 {
        return -5;
    }
    0
}