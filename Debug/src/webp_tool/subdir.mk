# 定义 webp_converter.o 的路径（相对于主 Makefile 所在的 Debug 目录）
src_webp_tool_webp_converter_O := src/webp_tool/webp_converter.o

# 定义 compress_40_quality_score.o 的路径
src_compress_40_quality_score_O := src/webp_tool/compress_40_quality_score.o

# 定义 webp_reconstructor_4.o 的路径
src_webp_tool_webp_reconstructor_4_O := src/webp_tool/webp_reconstructor_4.o

src_webp_tool_webp_reconstructor_40_O := src/webp_tool/webp_reconstructor_40.o

# 将目标文件添加到全局 OBJS 变量，确保主程序链接时包含它们
OBJS += $(src_webp_tool_webp_converter_O)
OBJS += $(src_compress_40_quality_score_O)
OBJS += $(src_webp_tool_webp_reconstructor_4_O)
OBJS += $(src_webp_tool_webp_reconstructor_40_O)

# 编译规则：将 webp_converter.cpp 编译为 .o 文件
$(src_webp_tool_webp_converter_O): src/webp_tool/webp_converter.cpp \
    src/webp_tool/webp_converter.h \
    makefile $(OPTIONAL_TOOL_DEPS)
	@echo 'Building file: $<'  # $< 表示依赖列表中的第一个文件（.cpp 源文件）
	@echo 'Invoking: Cross G++ Compiler'
	g++ $(CXXFLAGS) -MMD -MP -I/home/user/yexiang/libwebp/include \
		-c -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '

# 编译规则：将 compress_40_quality_score.cpp 编译为 .o 文件
$(src_compress_40_quality_score_O): src/webp_tool/compress_40_quality_score.cpp \
    src/webp_tool/compress_40_quality_score.h \
    makefile $(OPTIONAL_TOOL_DEPS)
	@echo 'Building file: $<'
	@echo 'Invoking: Cross G++ Compiler'
	g++ $(CXXFLAGS) -MMD -MP -I/home/user/yexiang/libwebp/include \
		-c -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '

# 编译规则：将 webp_reconstructor_4.cpp 编译为 .o 文件
$(src_webp_tool_webp_reconstructor_4_O): src/webp_tool/webp_reconstructor_4.cpp \
    src/webp_tool/webp_reconstructor_4.h \
    makefile $(OPTIONAL_TOOL_DEPS)
	@echo 'Building file: $<'
	@echo 'Invoking: Cross G++ Compiler'
	g++ $(CXXFLAGS) -MMD -MP -I/home/user/yexiang/libwebp/include \
		-c -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


$(src_webp_tool_webp_reconstructor_40_O): src/webp_tool/webp_reconstructor_40.cpp \
    src/webp_tool/webp_reconstructor_40.h \
    makefile $(OPTIONAL_TOOL_DEPS)
	@echo 'Building file: $<'  # $< 表示依赖列表中的第一个文件（.cpp 源文件）
	@echo 'Invoking: Cross G++ Compiler'
	g++ $(CXXFLAGS) -MMD -MP -I/home/user/yexiang/libwebp/include \
		-c -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '

# 引用生成的 .d 依赖文件（确保头文件修改时自动重新编译）
-include $(src_webp_tool_webp_converter_O:.o=.d)
-include $(src_compress_40_quality_score_O:.o=.d)
-include $(src_webp_tool_webp_reconstructor_4_O:.o=.d)
-include $(src_webp_tool_webp_reconstructor_40_O:.o=.d)