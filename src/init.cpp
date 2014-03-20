#include <Rcpp.h>


extern "C" {

#include <R_ext/Rdynload.h>

    RcppExport SEXP messina_messinaC(SEXP xSEXP, SEXP clsSEXP, SEXP n_bootSEXP, SEXP n_trainSEXP, SEXP minsensSEXP, SEXP minspecSEXP, SEXP progressSEXP, SEXP silentSEXP);


    R_CallMethodDef callMethods[] = {
        {"messina_messinaC", (DL_FUNC) &messina_messinaC, 8},
        {NULL, NULL, 0}
    };

    void R_init_messina(DllInfo *info)
    {
       // printf("hello from c++!\n"); // uncomment this to test
       R_registerRoutines(info, NULL, callMethods, NULL, NULL);
    }

}
