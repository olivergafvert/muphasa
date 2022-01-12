#ifndef MPH_SIGNATURECOLUMN_INCLUDED
#define MPH_SIGNATURECOLUMN_INCLUDED

#include "column.h"
#include "utils.h"

class signature_t : public std::pair<grade_t, index_t>{
 public:
    signature_t(grade_t _grade, index_t _index) : std::pair<grade_t, index_t>(_grade, _index){ }
    grade_t& get_grade() { return first; }
    index_t& get_index() { return second; }
};

template <typename Column> class BasicSignatureColumn : public Column{
public:
    grade_t grade;
    index_t signature_index;
    int local = 0; // Used in the presentation computation to classify column as local neg/pos or global.
    int last_updated;
    
    BasicSignatureColumn(const BasicSignatureColumn& column) : Column(column){
        this->grade = column.grade;
        this->signature_index = column.signature_index;
        this->last_updated = column.last_updated;
    }
    
    // TODO: Eventually, we can initilize these objects from the input buffers
    
    BasicSignatureColumn(const grade_t& grade, index_t signature_index) : Column(){
        this-> grade = grade;
        this->signature_index = signature_index;
        this->last_updated = 0;
    }
    
    BasicSignatureColumn(const grade_t& grade, index_t signature_index, int n) : Column(n){
        this-> grade = grade;
        this->signature_index = signature_index;
        this->last_updated = 0;
    }
    
    BasicSignatureColumn(const grade_t& grade, index_t signature_index, Column& column) : Column(column){
        this-> grade = grade;
        this->signature_index = signature_index;
        this->last_updated = 0;
    }
    
    grade_t get_grade() const { return grade; }
    
    
    void plus(BasicSignatureColumn& other_column){
        Column::plus(other_column);
        for(size_t i=0; i<this->grade.size(); i++){
            if(this->grade[i] < other_column.grade[i]){
                this->grade[i] = other_column.grade[i];
            }
        }
        return;
    }
    
    void plus_popped(BasicSignatureColumn& other_column){
        Column::plus_popped(other_column);
        for(size_t i=0; i<this->grade.size(); i++){
            if(this->grade[i] < other_column.grade[i]){
                this->grade[i] = other_column.grade[i];
            }
        }
        return;
    }
    
    void set_last_updated(size_t& index){
        last_updated = index;
    }
};


template <typename Column, typename Syzygy> class CustomSignatureColumn : public BasicSignatureColumn<Column> {
    /** Represents a graded column with signature information. */
public:
    Syzygy syzygy; // Keeps track of operations performed on the column.
    
    CustomSignatureColumn(const CustomSignatureColumn& column) : BasicSignatureColumn<Column>(column){
        this->syzygy = column.syzygy;
    }
    
    // TODO: Eventually, we can initilize these objects from the input buffers
    
    CustomSignatureColumn(const grade_t& grade, index_t signature_index) : BasicSignatureColumn<Column>(grade, signature_index){
        this->syzygy = Syzygy();
    }
    
    CustomSignatureColumn(const grade_t& grade, index_t signature_index, int n) : BasicSignatureColumn<Column>(grade, signature_index){
        this->syzygy = Syzygy(n);
    }
    
    CustomSignatureColumn(const grade_t& grade, index_t signature_index, int n, int m) : BasicSignatureColumn<Column>(grade, signature_index, n){
        this->syzygy = Syzygy(m);
    }
    
    CustomSignatureColumn(const grade_t& grade, index_t signature_index, Column& column) : BasicSignatureColumn<Column>(grade, signature_index, column){
        this->syzygy = Syzygy();
    }
    
    CustomSignatureColumn(const grade_t& grade, index_t signature_index, int n, Column& column) : BasicSignatureColumn<Column>(grade, signature_index, column){
        this->syzygy = Syzygy(n);
    }
    
    
    void plus(CustomSignatureColumn& other_column, bool should_add_syzygy=true){
        BasicSignatureColumn<Column>::plus(other_column);
        if(should_add_syzygy){
            this->syzygy.plus(other_column.syzygy);
        }
        return;
    }
    
    void plus_popped(CustomSignatureColumn& other_column){
        BasicSignatureColumn<Column>::plus_popped(other_column);
        this->syzygy.plus(other_column.syzygy);
        return;
    }
};


typedef VectorColumn SyzColumn;
typedef CustomSignatureColumn<VectorColumn, SyzColumn> SignatureColumn;

/*class SignatureColumnDense : public BasicSignatureColumn<DenseColumn> {
    // Represents a graded column with signature information.
public:
    DenseColumn syzygy; // Keeps track of operations performed on the column.
    
    SignatureColumnDense(const SignatureColumnDense& column) : BasicSignatureColumn<DenseColumn>(column){
        this->syzygy = column.syzygy;
    }
    
    // TODO: Eventually, we can initilize these objects from the input buffers
    
    SignatureColumnDense(const grade_t& grade, index_t signature_index) : BasicSignatureColumn<DenseColumn>(grade, signature_index){
        this->syzygy = DenseColumn();
    }
    
    SignatureColumnDense(const grade_t& grade, index_t signature_index, int n) : BasicSignatureColumn<DenseColumn>(grade, signature_index){
        this->syzygy = DenseColumn(n);
    }
    
    SignatureColumnDense(const grade_t& grade, index_t signature_index, DenseColumn& column) : BasicSignatureColumn<DenseColumn>(grade, signature_index, column){
        this->syzygy = DenseColumn();
    }
    
    SignatureColumnDense(const grade_t& grade, index_t signature_index, int n, DenseColumn& column) : BasicSignatureColumn<DenseColumn>(grade, signature_index, column){
        this->syzygy = DenseColumn(n);
    }
    
    void plus(SignatureColumnDense& other_column);
};


inline void SignatureColumnDense::plus(SignatureColumnDense& other_column)
{
    BasicSignatureColumn<DenseColumn>::plus(other_column);
    this->syzygy.plus(other_column.syzygy);
    return;
}
*/

// Tests
/*void column_addition_test(){
    spdlog::info("Running column addition test...");
    SignatureColumn column1(grade_t({1, 1, 0}), 0), column2(grade_t({0, 2, 1}), 1);
    
}*/

#endif // MPH_SIGNATURECOLUMN_INCLUDED
