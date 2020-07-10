
        #include <TMB.hpp>
        template<class Type>
        Type objective_function<Type>::operator() ()
        {
        PARAMETER_VECTOR(log_a);
        PARAMETER_VECTOR(log_b);
        PARAMETER_VECTOR(log_sigma);
        
        DATA_MATRIX(local_d);
        DATA_VECTOR(global_d);
        DATA_INTEGER(n_grid);
        DATA_VECTOR(grid);
        DATA_INTEGER(nyr);

        //int n = local_d.size();
        
        vector<Type> a = exp(log_a);
        vector<Type> b = exp(log_b);
        vector<Type> sigma = exp(log_sigma);
        vector<Type> pred_d(n_grid);
        //vector<Type> resid(n_grid);
        vector<Type> resid_quick(n_grid);
        matrix<Type> resid(n_grid,nyr);
        
        Type nll = 0.0;
        Type zero = 0.0;
        using namespace density; //package for various multivariate Gaussian distribs
        
      for (int j=0; j<n_grid; j++){
        vector<Type> quick_local_d = local_d.row(j);

        pred_d = a(j)*pow(global_d, b(j));
        resid_quick = quick_local_d-pred_d;
        resid.row(j) = resid_quick;
        
        nll -= dnorm(resid_quick, zero, sigma(j), true).sum();
      }

        REPORT(log_a);
        REPORT(log_b);
        REPORT(log_sigma);
        REPORT(resid);
        ADREPORT(log_b);
        return(nll);
        }
        
