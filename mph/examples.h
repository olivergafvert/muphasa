//
//  examples.h
//  mph
//
//  Created by Oliver on 2020-10-26.
//  Copyright © 2020 Oliver. All rights reserved.
//

#ifndef examples_h
#define examples_h
#include <iostream>
#include <chrono>
#include <random>
#include "utils.h"
#include "grade.h"
#include "signatureColumn.h"
#include "IO.h"

grade_t get_unique_grade(int p, std::vector<hash_set<index_t>>& visited_grades){
    grade_t grade;
    for(size_t i=0; i<p; i++){
        while(true){
            int val = rand()%10000;
            if(visited_grades[i].find(val) == visited_grades[i].end()){
                grade.push_back(val);
                visited_grades[i].insert(val);
                break;
            }
        }
    }
    return grade;
}

grade_t get_random_grade(int p, std::uniform_int_distribution<int>& dist, std::mt19937& seq){
    grade_t grade;
    for(size_t i=0; i<p; i++){
        int val = dist(seq);
        grade.push_back(val);
    }
    return grade;
}

void get_random_column(int n, SignatureColumn& column, double density){
    grade_t grade;
    double lower_bound = 0;
    double upper_bound = 1;
    std::uniform_real_distribution<double> unif(lower_bound,upper_bound);
    std::default_random_engine re;
    for(size_t i=0; i<n; i++){
        if(unif(re)< density){
            column.push(column_entry_t(1, i));
        }
    }
    return;
}

std::pair<Matrix, std::vector<grade_t>> random_boundary_matrix(int n, int m, int dim, int n_parameters, uint32_t seed=-1){
    /*
     Returns a random subset of the n-simplex boundary matrix at dimension 'dim' with random grades on the (dim-1)-simplices.
     */
    Matrix boundary;
    
    std::random_device rd;
    if(seed==-1){
        seed = rd() ^ (
                       (std::mt19937::result_type)std::chrono::duration_cast<std::chrono::seconds>(
                                                                                                   std::chrono::system_clock::now().time_since_epoch()
                                                                                                   ).count() +
                       (std::mt19937::result_type)
                       std::chrono::duration_cast<std::chrono::microseconds>(
                                                                             std::chrono::high_resolution_clock::now().time_since_epoch()
                                                                             ).count() );
    }
    
    std::mt19937 gen(seed);
    std::uniform_int_distribution<int> uni(0,m-1);
    
    std::vector<grade_t> row_grades;
    for(size_t i=0; i<m; i++){
        row_grades.push_back(get_random_grade(n_parameters, uni, gen));
    }
    std::vector<grade_t> point_grades;
    for(size_t i=0; i<m; i++){
        grade_t column_grade = get_random_grade(n_parameters, uni, gen);
        point_grades.push_back(row_grades[i].join(column_grade));
    }
    for(size_t i=0; i<n; i++){
        hash_set<int> visited;
        int index = uni(gen);
        visited.insert(index);
        SignatureColumn column(point_grades[index], i);
        column.push(column_entry_t(1, index));
        grade_t row_grade = row_grades[index];
        for(size_t j=0; j<dim-1; j++){
            index = uni(gen);
            while(visited.find(index) != visited.end()){
                index = uni(gen);
            }
            visited.insert(index);
            column.push(column_entry_t(1, index));
            row_grade = row_grade.join(row_grades[index]);
            //column.grade = column.grade.join(point_grades[index]);
        }
        column.grade = row_grade;
        boundary.push_back(column);
    }
    return std::pair<Matrix, std::vector<grade_t>>(boundary, row_grades);
}


std::vector<input_t> get_random_next_step(std::vector<input_t> point, std::normal_distribution<double>& dist, std::mt19937& seq){
    /* Assume last coordinate is time */
    std::vector<input_t> point_;
    for(size_t i=0; i<point.size()-1; i++){
        point_.push_back(point[i] + dist(seq));
    }
    point_.push_back(point[point.size()-1] + 1);
    return point_;
}

std::vector<std::vector<input_t>> trajectory(size_t n_samples, size_t dim, uint32_t seed=-1){
    std::random_device rd;
    if(seed==-1){
        seed = rd() ^ (
                                             (std::mt19937::result_type)std::chrono::duration_cast<std::chrono::seconds>(
                                                                                                                         std::chrono::system_clock::now().time_since_epoch()
                                                                                              ).count() +
                                             (std::mt19937::result_type)
                                             std::chrono::duration_cast<std::chrono::microseconds>(
                                                                                                   std::chrono::high_resolution_clock::now().time_since_epoch()
                                                                                                   ).count() );
    }
    
    std::mt19937 gen(seed);
    std::normal_distribution<double> d{0,1};
    std::uniform_real_distribution<double> uniform(0.0,1.0);
    
    std::vector<std::vector<input_t>> points;
    std::vector<input_t> point;
    for(size_t i=0; i<dim+1; i++){
        point.push_back(uniform(gen));
    }
    for(size_t i=0; i<n_samples; i++){
        point = get_random_next_step(point, d, gen);
        points.push_back(point);
    }
    return points;
}

std::vector<std::vector<input_t>> trajectory_box(size_t n_samples, size_t dim, uint32_t seed=-1){
    std::random_device rd;
    if(seed==-1){
        seed = rd() ^ (
                       (std::mt19937::result_type)std::chrono::duration_cast<std::chrono::seconds>(
                                                                                                   std::chrono::system_clock::now().time_since_epoch()
                                                                                                   ).count() +
                       (std::mt19937::result_type)
                       std::chrono::duration_cast<std::chrono::microseconds>(
                                                                             std::chrono::high_resolution_clock::now().time_since_epoch()
                                                                             ).count() );
    }
    
    std::mt19937 gen(seed);
    std::normal_distribution<double> d{0,0.2};
    std::uniform_real_distribution<double> uniform(0.0, 1.0);
    
    std::vector<std::vector<input_t>> points;
    std::vector<input_t> point;
    for(size_t i=0; i<dim+1; i++){
        point.push_back(uniform(gen));
    }
    for(size_t i=0; i<n_samples; i++){
        point = get_random_next_step(point, d, gen);
        for(size_t j=0; j<point.size(); j++){
            if(point[j] < 0){
                point[j] = -point[j];
            }else if(point[j] > 1){
                point[j] = 2-point[j];
            }
        }
        points.push_back(point);
    }
    return points;
}


std::vector<std::vector<input_t>> time_varying_point_cloud(size_t n_points, size_t n_timesteps, size_t dim, uint32_t seed=-1){
    std::vector<std::vector<input_t>> points;
    for(size_t i=0; i<n_points; i++){
        std::vector<std::vector<input_t>> t;
        if(seed==-1){
            t = trajectory_box(n_timesteps, dim);
        }else{
            t = trajectory_box(n_timesteps, dim, seed+(uint32_t)i);
        }
        for(auto& p:t){
            points.push_back(p);
        }
    }
    return points;
}

std::vector<std::vector<std::vector<input_t>>> time_varying_point_cloud_trajectories(size_t n_points, size_t n_timesteps, size_t dim, uint32_t seed=-1){
    std::vector<std::vector<std::vector<input_t>>> trajectories;
    for(size_t i=0; i<n_points; i++){
        std::vector<std::vector<input_t>> t;
        if(seed==-1){
            t = trajectory_box(n_timesteps, dim);
        }else{
            t = trajectory_box(n_timesteps, dim, seed+(uint32_t)i);
        }
        trajectories.push_back(t);
    }
    return trajectories;
}

std::pair<Matrix, Matrix> get_random_boundary_matrix(size_t n_points, size_t n_dim, size_t n_parameters, uint32_t seed=-1){
    std::random_device rd;
    if(seed==-1){
        seed = rd() ^ (
                       (std::mt19937::result_type)std::chrono::duration_cast<std::chrono::seconds>(
                                                                                                   std::chrono::system_clock::now().time_since_epoch()
                                                                                                   ).count() +
                       (std::mt19937::result_type)
                       std::chrono::duration_cast<std::chrono::microseconds>(
                                                                             std::chrono::high_resolution_clock::now().time_since_epoch()
                                                                             ).count() );
    }
    
    std::mt19937 gen(seed);
    std::uniform_int_distribution<int> uni(0,100);
    
    grade_t max_grade;
    for(size_t i=0; i<n_parameters; i++){
        max_grade.push_back(100);
    }
    
    
    std::vector<std::vector<input_t>> points;
    for(size_t i=0; i<n_points; i++){
        points.push_back(std::vector<input_t>{1});
    }
    std::vector<grade_t> grades;
    for(size_t i=0; i<n_points; i++){
        grades.push_back(get_random_grade(n_parameters, uni, gen));
    }
    return compute_boundary_matrices_grades(points, grades, max_grade, n_dim);
}

#endif /* examples_h */
