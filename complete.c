/*
 * complete_fast_bigint.c  –  O(n^3) DP with GMP arbitrary-precision integers.
 *
 * Compile:  gcc -O2 -o complete_fast_bigint complete_fast_bigint.c -lgmp
 *
 * DP recurrence (see complete_fast.c for full derivation):
 *
 *   dp[r][k][p] = number of ways to complete a bracket sequence that has
 *                 r pairs still to open, k intervals currently open, and
 *                 crossings-so-far ≡ p (mod 2), such that the TOTAL number
 *                 of crossings is even.
 *
 *   Base case:    dp[0][0][0] = 1,  dp[0][0][1] = 0
 *
 *   Transitions:
 *     Open a new pair  (r > 0):   += r * dp[r-1][k+1][p]
 *     Close (even-position slot):  += ceil(k/2)  * dp[r][k-1][p]
 *     Close (odd-position slot):   += floor(k/2) * dp[r][k-1][p^1]
 *
 * tot_num(n) = (2n)! / 2^n   (total constrained permutations)
 * bias       = 2*(good - tot_num/2) / tot_num
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <gmp.h>

#define MAXN 200   /* supports graph sizes up to 200; increase freely */

/* ── Memo table ──────────────────────────────────────────────────────────── */

static mpz_t memo[MAXN+1][MAXN+1][2];
static int   initialised[MAXN+1][MAXN+1][2];
static int   table_ready = 0;

static void init_table(void) {
    if (table_ready) return;
    for (int r = 0; r <= MAXN; r++)
        for (int k = 0; k <= MAXN; k++)
            for (int p = 0; p < 2; p++) {
                mpz_init(memo[r][k][p]);
                initialised[r][k][p] = 0;
            }
    table_ready = 1;
}

static void clear_table(int n) {
    for (int r = 0; r <= n; r++)
        for (int k = 0; k <= n; k++)
            for (int p = 0; p < 2; p++)
                initialised[r][k][p] = 0;
}

/* ── DP ──────────────────────────────────────────────────────────────────── */

static void dp(mpz_t result, int r, int k, int p) {
    if (r == 0 && k == 0) {
        mpz_set_ui(result, p == 0 ? 1 : 0);
        return;
    }

    if (initialised[r][k][p]) {
        mpz_set(result, memo[r][k][p]);
        return;
    }

    mpz_t tmp;
    mpz_init(tmp);
    mpz_set_ui(memo[r][k][p], 0);

    /* Open a new interval */
    if (r > 0) {
        dp(tmp, r - 1, k + 1, p);
        mpz_addmul_ui(memo[r][k][p], tmp, (unsigned long)r);
    }

    /* Close an open interval */
    if (k > 0) {
        unsigned long even_choices = (unsigned long)(k + 1) / 2;
        unsigned long odd_choices  = (unsigned long) k      / 2;

        if (even_choices > 0) {
            dp(tmp, r, k - 1, p);
            mpz_addmul_ui(memo[r][k][p], tmp, even_choices);
        }
        if (odd_choices > 0) {
            dp(tmp, r, k - 1, p ^ 1);
            mpz_addmul_ui(memo[r][k][p], tmp, odd_choices);
        }
    }

    mpz_clear(tmp);
    initialised[r][k][p] = 1;
    mpz_set(result, memo[r][k][p]);
}

/* ── tot_num = (2n)! / 2^n ───────────────────────────────────────────────── */

static void compute_tot_num(mpz_t out, int n) {
    mpz_fac_ui(out, (unsigned long)(2 * n));
    mpz_tdiv_q_2exp(out, out, (mp_bitcnt_t)n);
}

/* ── Main ────────────────────────────────────────────────────────────────── */

int main(void) {
    int max_graph_size = 50;

    init_table();

    mpz_t good, total, two_good, bias_num;
    mpf_t bias_f, total_f, num_f;

    mpz_inits(good, total, two_good, bias_num, NULL);
    mpf_set_default_prec(256);
    mpf_inits(bias_f, total_f, num_f, NULL);

    for (int size = 2; size <= max_graph_size; size++) {
        clear_table(size);

        dp(good, size, 0, 0);
        compute_tot_num(total, size);
        gmp_printf("Graph size: %d\n", size);
        gmp_printf("%Zd %Zd\n", good, total);

        /* bias = (2*good - total) / total */
        mpz_mul_ui(two_good, good, 2);
        mpz_sub(bias_num, two_good, total);

        mpf_set_z(num_f,   bias_num);
        mpf_set_z(total_f, total);
        mpf_div(bias_f, total_f, num_f);

        gmp_printf("Inverse of expected value: %.15Ff\n", bias_f);
    }

    mpz_clears(good, total, two_good, bias_num, NULL);
    mpf_clears(bias_f, total_f, num_f, NULL);
    return 0;
}