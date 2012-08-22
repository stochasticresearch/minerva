#include <stdlib.h>


typedef struct mine_problem
{
  int n;
  double *x;
  double *y;
} mine_problem;

typedef struct mine_parameter
{
  double alpha;
  double c;
} mine_parameter;

typedef struct mine_score
{
  int m;
  int *p;
  double **I;
} mine_score;

mine_score *mine_compute_score(mine_problem *prob, mine_parameter *param);
char *check_parameter(mine_parameter *param);
double mic(mine_score *score);
double mas(mine_score *score);
double mev(mine_score *score);
double mcn(mine_score *score);
void mine_free_score(mine_score **score);
