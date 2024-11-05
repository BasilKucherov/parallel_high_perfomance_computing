namespace MathUtils {
void AdamarsMult(int N, const double *X, const double *Y, double *Res);
void SpMV(int N, const int *IA, const int *JA, const double *A, const double *X,
          double *Res);
void DotProduct(int N, const double *X, const double *Y, double &Res);
void Axpy(int N, double a, const double *X, const double *Y, double *Res);

void Solve(int N, const int *IA, const int *JA, const double *A,
           const double *B, double Eps, int MaxIter, double *&X, int &Iter,
           double &Res, bool Debug);
} // namespace MathUtils
