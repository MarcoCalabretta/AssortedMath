const double error = 0.00001;
int doub_eq(double a, double b){
    double difference = a - b;
    return difference < error && -difference < error;
}

void print_matrix(int rows, int cols, double matrix[rows][cols]){
    for(int i = 0; i < rows; i++){
        for(int j = 0; j < cols; j++){
            printf("%.03f, ", matrix[i][j]);
        }
        printf("\n");
    }
    printf("\n");
}

void row_swap(int i, int j, int rows, int cols, double matrix[rows][cols]){
    double temp;
    for(int count = 0; count < cols; count++){
        temp = matrix[j][count];
        matrix[j][count] = matrix[i][count];
        matrix[i][count] = temp;
    }
}
void row_add(int actor, double scalar, int acted, int rows, int cols, double matrix[rows][cols]){
    for(int count = 0; count < cols; count++){
        matrix[acted][count] += matrix[actor][count]*scalar;
    }
}
void row_scale(int i, double scalar, int rows, int cols, double matrix[rows][cols]){
    for(int count = 0; count < cols; count++){
        matrix[i][count] *= scalar;
    }
}

// checks if a matrix is in row reduced echelon form (rref)
// uses 4 conditions:
// 1. All rows containing a non-zero entry are above rows which only contains zeros.
// 2. The first non-zero entry in each non-zero row is 1, called a leading one.
// 3. The leading one in each non-zero row is to the right of the leading one in any
// row above it.
// 4. A leading one is the only non-zero entry in its column.
// returns 0 if matrix is in rref, otherwise it returns the number of condition that was violated
int rref_check(int rows, int cols, double matrix[rows][cols]){
    int zero_row = 0;
    int leading_column = -1;

    for(int i = 0; i < rows; i++){
        for(int j = 0; j < cols; j++){
            if(!doub_eq(matrix[i][j], 0)){

                //condition 1 checks if all non-zero rows are above all zero rows
                //if there was already a zero row and we find a non-zero value it's not in rref
                if(zero_row) return 1;

                //condition 2 checks if the first leading entry is a 1
                if(!doub_eq(matrix[i][j], 1)) return 2;

                // condition 3 checks if the leading 1 is to the right of other leading ones
                if(j <= leading_column){
                    return 3;
                } else{
                    leading_column = j;
                }
                //condition 4 checks if the leading 1 is the only entry in its column
                for(int k = 0; k < cols; k++){
                    if(k != i && !doub_eq(matrix[k][j], 0)){
                        return 4;
                    }
                }

                break;
            }
            //if it gets to the end of the row only finding zeros it sets zero_row to 1;
            if(j == cols - 1){
                zero_row = 1;
            }
        }
    }
    return 0;
}

void gaussian_elimination(int rows, int cols, double matrix[rows][cols]){
    int rref = -1;
    while(rref != 0){
        rref = rref_check(rows, cols, matrix);
        if(rref == 0){
            break;
        } else if(rref == 1){
            // find first zero row
            int zero_row = -1;
            for(int i = 0; i < rows; i++){
                for(int j = 0; j < cols; j++){
                    if(!doub_eq(matrix[i][j], 0)) break;
                    if(j == cols - 1) zero_row = i;
                }
                if(zero_row != -1) break;
            }

            // find the last non-zero row
            int non_zero_row = -1;
            for(int i = rows - 1; i >= 0; i--){
                for(int j = 0; j < cols; j++){
                    if(!doub_eq(matrix[i][j], 0)){
                        non_zero_row = i;
                        break;
                    }
                }
                if(non_zero_row != -1) break;
            }
            // swap them
            row_swap(zero_row, non_zero_row, rows, cols, matrix);

        } else if(rref == 2){
            // find first leading value that's not a 1
            double scalar = 1;
            int scaled_row = -1;
            for(int i = 0; i < rows; i++){
                for(int j = 0; j < cols; j++){
                    if(!doub_eq(matrix[i][j], 0)){
                        if(!doub_eq(matrix[i][j], 1)){
                            scalar = 1/matrix[i][j];
                            scaled_row = i;
                        }
                        break;
                    }
                }
                if(scalar != 1) break;
            }

            // make it a 1
            row_scale(scaled_row, scalar, rows, cols, matrix);

        } else if(rref == 3){
            // find first row that has a leading 1 to the left
            int leading_column = -1;
            int out_of_order_row = -1;
            for(int i = 0; i < rows; i++){
                for(int j = 0; j < cols; j++){
                    if(j <= leading_column){
                        out_of_order_row = i;
                        break;
                    } else{
                        leading_column = j;
                    }
                }
                if(out_of_order_row != -1) break;
            }

            // switch it with the row on top
            row_swap(out_of_order_row, out_of_order_row + 1, rows, cols, matrix);

        } else if(rref == 4){
            // find row with leading 1 that's not alone in its column
            int off_column = -1;
            int off_row = -1;
            // find non-zero column values
            double column_values[cols];
            for(int i = 0; i < rows; i++){
                for(int j = 0; j < cols; j++){
                    if(!doub_eq(matrix[i][j], 0)){
                        for(int k = 0; k < cols; k++){
                            column_values[k] = matrix[k][j];
                            if(k != i && !doub_eq(matrix[k][j], 0)){
                                off_column = j;
                                off_row = i;
                            }
                        }
                        break;
                    }
                }
                if(off_column != -1) break;
            }
            // for each non-zero value, add copies of the original row to make it a zero value
            for(int i = 0; i< cols; i++){
                if(i != off_row && !doub_eq(column_values[i], 0)){
                    row_add(off_row, -column_values[i]/column_values[off_row], i, rows, cols, matrix);
                }
            }
        }
    }
}

void main(){
    int rows = 9, cols = 15;
    double matrix[9][15] = {
    {30, 51, 54, 26, 48, 67, 97, 13, 81, 90, 88, 58, 14, 14, 80},
    {83, 81, 15, 2, 41, 4, 83, 79, 69, 66, 63, 35, 41, 18, 87},
    {47, 67, 82, 51, 11, 30, 41, 90, 57, 47, 33, 33, 34, 17, 87},
    {32, 19, 80, 64, 58, 70, 99, 73, 15, 42, 72, 51, 82, 56, 55},
    {70, 13, 64, 17, 42, 18, 23, 68, 20, 57, 94, 33, 86, 6, 88},
    {99, 72, 66, 82, 4, 99, 14, 57, 10, 83, 62, 98, 55, 59, 16},
    {14, 95, 38, 83, 35, 54, 81, 49, 23, 91, 5, 30, 18, 100, 26},
    {60, 15, 25, 2, 60, 79, 38, 17, 56, 94, 41, 46, 75, 25, 62},
    {76, 6, 38, 2, 56, 47, 9, 24, 62, 80, 34, 17, 1, 40, 35},
    };
    gaussian_elimination(rows, cols, matrix);
    print_matrix(rows, cols, matrix);
}
