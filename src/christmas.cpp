#include <iostream>
#include <unistd.h>   // Linux 延时函数 usleep
#include <cstdlib>    // 随机数
#include <ctime>      // 随机数种子
#include <sys/ioctl.h> // 获取终端宽度（可选，优化居中）
#include <unistd.h>
using namespace std;

// Linux 控制台文字颜色设置（ANSI 转义序列）
// color: 1-7（基础色）/90-97（亮色），参考 ANSI 颜色码
void setColor(int color) {
    // 基础色（30-37）：黑、红、绿、黄、蓝、紫、青、白
    // 亮色（90-97）：亮黑、亮红、亮绿、亮黄、亮蓝、亮紫、亮青、亮白
    int ansiColor = 30 + (color % 8);  // 基础色范围
    if (color >= 8) {
        ansiColor = 90 + (color % 8); // 亮色范围
    }
    cout << "\033[" << ansiColor << "m"; // ANSI 颜色控制码
}

// 恢复终端默认颜色
void resetColor() {
    cout << "\033[0m";
}

// 获取终端宽度（可选，优化圣诞树居中）
int getTerminalWidth() {
    struct winsize w;
    ioctl(STDOUT_FILENO, TIOCGWINSZ, &w);
    return w.ws_col;
}

int main() {
    srand(time(0)); // 初始化随机数种子（避免每次彩灯颜色相同）
    int height = 10; // 圣诞树高度，可调整（建议5-15）
    int termWidth = getTerminalWidth(); // 获取终端宽度，优化居中

    // 圣诞树主体（三角形部分）
    for (int i = 0; i < height; i++) {
        int starCount = 2 * i + 1;
        // 计算前置空格（基于终端宽度居中，更适配不同终端）
        int spaceCount = (termWidth - starCount) / 2;
        // 打印前置空格（居中效果）
        for (int j = 0; j < spaceCount; j++) {
            cout << " ";
        }
        // 打印圣诞树分支（带随机彩灯）
        for (int k = 0; k < starCount; k++) {
            int randColor = rand() % 16; // 0-15 随机颜色（覆盖基础色+亮色）
            setColor(randColor);
            cout << "*";
        }
        resetColor(); // 临时重置，避免颜色溢出到换行
        cout << endl;
        usleep(100000); // 每一层延时100ms（Linux 中 usleep 单位是微秒，100ms=100000μs）
    }

    // 圣诞树树干（矩形部分）
    int trunkWidth = 3; // 树干宽度
    for (int i = 0; i < 3; i++) {  // 树干高度3层，可调整
        // 树干居中
        int spaceCount = (termWidth - trunkWidth) / 2;
        for (int j = 0; j < spaceCount; j++) {
            cout << " ";
        }
        setColor(33); // 树干固定棕色（ANSI 33=黄色/棕色，93=亮黄）
        cout << "###" << endl;  // 树干宽度3个#，可调整
        resetColor();
        usleep(100000);
    }

    // 底部装饰文字
    string msg = "Merry Christmas!";
    int msgSpace = (termWidth - msg.length()) / 2; // 文字居中
    setColor(93); // 亮黄色文字（ANSI 93）
    for (int i = 0; i < msgSpace; i++) {
        cout << " ";
    }
    cout << msg << endl;

    resetColor(); // 恢复终端默认颜色
    return 0;
}