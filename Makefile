OBJS =  main.o cal_Fai.o cal_Jz.o cal_Jz_coef.o cal_RHS_vector.o sum.o cal_Ez.o cal_Exp.o cal_Eyp.o cal_Ezp.o cal_Hxp.o cal_Hyp.o cal_Hzp.o cal_Hzxp.o cal_Hzyp.o cal_Fxp.o cal_Fyp.o cal_Fzxp.o cal_Fzyp.o cal_Jmxp.o cal_Jmyp.o cal_Jmzp.o get_node_number.o compose_coef_matrix.o cal_Nyu.o cal_P.o memory_allocate3d.o memory_allocate3d_int.o memory_allocate2d.o memory_allocate1d.o delete3d.o delete3d_int.o delete2d.o delete1d.o PML_cond.o

main: $(OBJS)
	g++ $(OBJS) -o main

%.o: %.cpp WLP2D.h
	g++ -c $< -Wall
